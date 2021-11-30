#!/usr/bin/env python3

import stepwise
import stepwise_mol_bio
import freezerbox
import pytest
import parametrize_from_file

from voluptuous import Schema, Invalid, Coerce, And, Or, Optional
from unittest.mock import MagicMock
from contextlib import nullcontext
from textwrap import indent

class do_with:
    # Parent class for `eval_with` and `exec_with`.  Not meant to be used 
    # directly.

    def __init__(self, globals=None, **kw_globals):
        self.globals = globals or {}
        self.globals.update(kw_globals)

    def use(self, **kw_globals):
        self.globals.update(kw_globals)
        return self

    def all(self, module):
        try:
            keys = module.__all__
        except AttributeError:
            keys = module.__dict__

        self.globals.update({
                k: module.__dict__[k]
                for k in keys
        })
        return self

class eval_with(do_with):
    # Evaluate snippets of python code, with powerful control over the names 
    # available to that snippet.
    # 
    # Examples:
    #
    # With no arguments, this just creates a callable that will evaluate 
    # snippets of python code:
    #
    #     >>> eval_with()('1 + 1')
    #     2
    #
    # With arguments, you can add global variables to the namespace where the 
    # snippet is evaluated.  This is most often use to "import" modules:
    #
    #    >>> import math
    #    >>> eval_math = eval_with(math=math)
    #    >>> eval_math('math.sqrt(2)')
    #    1.4142135623730951
    #
    # You can use the fluent-style `use()` method to add global variables after 
    # the object has been created.  This is mostly useful for building on 
    # existing evaluators:
    #
    #    >>> import numpy as np
    #    >>> eval_np = eval_math.use(np=np)
    #    >>> eval_np('np.arange(3)')
    #    array([0, 1, 2])
    #
    # The `all()` method provides a shorthand way to import all names from a 
    # module:
    #
    #    >>> eval_math = eval_with().all(math)
    #    >>> eval_math('sqrt(2)')
    #    1.4142135623730951
    #
    # Details:
    #
    # The input to evaluate can be a list, dict, or string.  Strings are 
    # directly evaluated.  List items and dict values are recursively 
    # evaluated.  Dict keys are not processed.  This allows you to switch 
    # freely between nested test and python syntax, depending on which makes 
    # most sense for each particular input.
    #
    # This method is meant to be used with voluptuous, as it raises 
    # `voluptuous.Invalid` in the event that the given input can't be 
    # evaluated.

    def __call__(self, src):
        try:
            if isinstance(src, list):
                return [self(x) for x in src]
            elif isinstance(src, dict):
                return {k: self(v) for k, v in src.items()}
            else:
                return eval(src, self.globals)

        except Exception as err:
            raise Invalid(str(err)) from err

class exec_with(do_with):
    # This has the same API as `eval_with`, except the given snippet is exec'd 
    # instead of eval'd.  You must specify a name, and the value of the 
    # variable with that name in the exec'd namespace is what will be 
    # returned.
    # 
    # Example:
    # 
    #     >>> f = exec_with('x')
    #     >>> f('x = 1')
    #     1

    def __init__(self, get, globals=None, **kw_globals):
        super().__init__(globals, **kw_globals)
        self.get = get

    def get(self, get):
        self.get = get

    def __call__(self, src):
        globals = self.globals.copy()

        try:
            exec(src, globals)
        except Exception as err:
            raise Invalid(str(err)) from err

        if callable(self.get):
            return self.get(globals)
        else:
            return globals[self.get]

def error_or(expected):
    schema = {}

    # Either specify an error or an expected value, not both.
    # KBK: This doesn't work for some reason.
    #schema[Or('error', *expected, only_one=True)] = object

    schema[Optional('error', default='none')] = error

    def check(schema):
        # This function basically just reimplements `Or`, but with a better 
        # error message.  The error message for `Or` is confusing because it 
        # only mentions the first option and not the second.  That's especially 
        # bad in this case, where the first option is basically an internal 
        # implementation detail.  This function changes the error message so 
        # that only the second option is mentioned, which is the least 
        # confusing in this case.
        # 
        # It's worth noting that the real problem is that voluptuous checks 
        # default values.  I can't imagine a case where this is actually 
        # useful. I might think about making a PR for that.

        def do_check(x):
            if isinstance(x, MagicMock):
                return x
            return Schema(schema)(x)
        return do_check

    schema.update({
        Optional(k, default=MagicMock()): check(v)
        for k, v in expected.items()
    })
    return schema

# Something to think about: I'd like to put a version of this function in the 
# `parametrize_from_file` package.  I need a general way to specify the local 
# variables, though.  And to `eval()` the exception type...

def error(x):
    if x == 'none':
        return nullcontext()

    err_eval = eval_with(
            ConfigError=stepwise_mol_bio.ConfigError,
            UsageError=stepwise_mol_bio.UsageError,
            QueryError=freezerbox.QueryError,
    )
    if isinstance(x, str):
        err_type = err_eval(x)
        err_messages = []
    else:
        err_type = err_eval(x['type'])
        err_messages = x.get('message', [])
        if not isinstance(err_messages, list):
            err_messages = [err_messages]

    # Normally I'd use `@contextmanager` to make a context manager like this, 
    # but generator-based context managers cannot be reused.  This is a problem 
    # for tests, because if a test using this context manager is parametrized, 
    # the same context manager instance will need to be reused multiple times.  
    # The only way to support this is to implement the context manager from 
    # scratch.

    class expect_error:

        def __repr__(self):
            return f'<expect_error type={err_type!r} messages={err_messages!r}>'

        def __enter__(self):
            self.raises = pytest.raises(err_type)
            self.err = self.raises.__enter__()

        def __exit__(self, *args):
            if self.raises.__exit__(*args):
                for msg in err_messages:
                    self.err.match(msg)
                return True

    return expect_error()

def empty_ok(x):
    return Or(x, And('', lambda y: type(x)()))

def eval_db(reagents):
    db = freezerbox.Database({'use': 'TEST_DB'})
    reagents = Schema(empty_ok({str: str}))(reagents)

    eval_local = eval_with(DB=db, **eval_swmb.globals)

    for tag, reagent in reagents.items():
        db[tag] = eval_local(reagent)

    return db

def exec_app(src):
    app = exec_with('app', stepwise=stepwise)\
            .all(stepwise)\
            .all(freezerbox)\
            .all(stepwise_mol_bio)(src)

    return app

eval_py = eval_with()
eval_swmb = eval_with().all(stepwise).all(freezerbox).all(stepwise_mol_bio)

app_expected = Schema({
    'app': exec_app,
    'expected': eval,
})
app_expected_reaction = Schema({
    'app': exec_app,
    'expected': str,
})
app_expected_protocol = Schema({
    'app': exec_app,
    'expected': [str],
})
app_expected_protocol_error = Schema({
    'app': exec_app,
    **error_or({
        'expected': [str],
    }),
})
app_expected_error = Schema({
    'app': exec_app,
    **error_or({
        'expected': eval,
    }),
})
db_expected = Schema({
    'db': eval_db,
    'expected': eval_swmb,
})
db_expected_protocol = Schema({
    'db': eval_db,
    'expected': [str],
})
db_tags_expected_protocol = Schema({
    'db': eval_db,
    'tags': [str],
    'expected': [str],
})
cmd_stdout_stderr = Schema({
    'cmd': str,
    Optional('stdout', default='^$'): str,
    Optional('stderr', default='^$'): str,
})

def match_protocol(app, expected, capture=nullcontext()):
    with capture:
        actual = app.protocol.format_text()

    prev = None
    print(actual.strip() + '\n')

    i = 0
    for x in expected:
        j = actual.find(x, i)

        if j == -1:
            print(f"Expected:\n  {x!r}")
            if prev:
                print(f"After:\n  {prev!r}")

            return False

        i = j + len(x)
        prev = x

    return True
