import stepwise
import stepwise_mol_bio
import freezerbox
import pytest
import parametrize_from_file

from voluptuous import Schema, Invalid, Coerce, And, Or, Optional
from parametrize_from_file.voluptuous import Namespace, empty_ok
from stepwise.testing import check_command, disable_capture
from freezerbox.stepwise import Make
from more_itertools import always_iterable
from contextlib import nullcontext

with_py = Namespace()
with_swmb = Namespace(
        'import stepwise',
        'from stepwise import *',
        'from freezerbox import *',
        'from stepwise_mol_bio import *',
)

def eval_db(reagents):
    db = freezerbox.Database({'use': 'TEST_DB'})
    reagents = Schema(empty_ok({str: str}))(reagents)

    with_local = Namespace(with_swmb, DB=db)

    for tag, reagent in reagents.items():
        db[tag] = with_local.eval(reagent)

    return db

def exec_app(src):
    return with_swmb.exec(src, get='app')

def match_protocol(app, expected, forbidden=[], *, capture=nullcontext()):
    with capture:
        actual = app.protocol.format_text()

    prev = None
    print(actual.strip() + '\n')

    i = 0
    for x in always_iterable(expected):
        j = actual.find(x, i)

        if j == -1:
            print(f"Expected:\n  {x!r}")
            if prev:
                print(f"After:\n  {prev!r}")

            return False

        i = j + len(x)
        prev = x

    for x in always_iterable(forbidden):
        if x in actual:
            print(f"Didn't expect:\n  {x!r}")
            return False

    return True

def noop(x):
    return x


def parametrize_from_file_factory(*args, **kwargs):
    """
    Return a factory function that, when invoked, will create a copy of the 
    decorated test function that will read test parameters from the data file 
    corresponding to script that the factory was invoked from.  In other words, 
    if the factory is invoked in `test_spam.py`, it will create a test function 
    that reads parameters from `test_spam.nt`.  The new test function is also 
    automatically added to the caller's global scope, so that pytest will 
    be able to find it.

    The purpose of making factories like these is to use the same test function 
    for multiple protocols.
    """
    import inspect, sys
    from pathlib import Path

    def decorator(f):

        def factory():
            frame = inspect.currentframe()
            try:
                ns = frame.f_back.f_globals
                path_py = Path(ns['__file__'])

                test_func = copy_func(f)
                test_func.__module__ = path_py.stem
                test_params = parametrize_from_file(*args, **kwargs)
                ns[test_func.__name__] = test_params(test_func)

            finally:
                del frame

        factory.__test__ = False
        return factory

    def copy_func(f, name=None):
        # https://stackoverflow.com/questions/6527633/how-can-i-make-a-deepcopy-of-a-function-in-python
        from types import FunctionType
        from copy import copy
        g = FunctionType(
                f.__code__,
                f.__globals__,
                name or f.__name__,
                f.__defaults__,
                f.__closure__,
        )
        g.__dict__ = copy(f.__dict__)
        return g

    return decorator

@parametrize_from_file_factory()
def test_reaction(app, expected):
    app = exec_app(app)
    assert app.reaction.format_text() == expected

@parametrize_from_file_factory(
        schema=Schema({
            'app': str,
            **with_swmb.error_or({
                'expected': Or(str, [str]),
                'forbidden': Or(str, [str]),
            }),
        }),
)
def test_protocol(app, expected, forbidden, error):
    app = exec_app(app)
    with error:
        assert match_protocol(app, expected, forbidden)

@parametrize_from_file_factory(
        schema=Schema({
            'cmd': str,
            Optional('stdout', default='^$'): str,
            Optional('stderr', default='^$'): str,
        }),
)
@pytest.mark.slow
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout, stderr=stderr)

@parametrize_from_file_factory(
        schema=Schema({
            'db': {str: str},
            Optional('tags', default=[]): [str],
            'expected': [str],
        }),
)
def test_freezerbox_make(db, tags, expected, disable_capture):
    db = eval_db(db)
    tags = tags or list(db.keys())
    app = Make(db, tags)
    assert match_protocol(app, expected, capture=disable_capture)

@parametrize_from_file_factory(
        schema=Schema({
            'db': {str: str},
            'expected': {str: {str: str}},
            Optional('errors', default={}): {str: {str: str}}
        }),
)
def test_freezerbox_attrs(db, expected, errors):
    db = eval_db(db)
    expected = with_swmb.eval(expected)

    for tag in expected:
        for attr in expected[tag]:
            assert getattr(db[tag], attr) == expected[tag][attr]

    for tag in errors:
        for attr in errors[tag]:
            with with_swmb.error(errors[tag][attr]):
                getattr(db[tag], attr)

