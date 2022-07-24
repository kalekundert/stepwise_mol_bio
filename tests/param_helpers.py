import stepwise
import stepwise_mol_bio
import freezerbox
import pytest
import parametrize_from_file

from voluptuous import Schema, Invalid, Coerce, And, Or, Optional
from parametrize_from_file import Namespace, error_or, defaults, cast
from stepwise_mol_bio import make_namespace_recursive
from stepwise.testing import check_command, disable_capture
from pytest import approx
from pytest_unordered import unordered
from freezerbox.stepwise import Make
from more_itertools import always_iterable
from contextlib import nullcontext
from difflib import Differ
from types import SimpleNamespace

with_py = Namespace()
with_math = Namespace('from math import *')
with_sw = Namespace(
        'import stepwise',
        'from stepwise import *',
)
with_swmb = Namespace(
        with_sw,
        'from freezerbox import *',
        'from stepwise_mol_bio import *',
)

Int = Coerce(int)
Float = Coerce(float)

def eval_db(reagents):
    db = freezerbox.Database({'use': 'TEST_DB'})
    reagents = Schema(empty_ok({str: str}))(reagents)

    with_local = Namespace(with_swmb, DB=db)

    for tag, reagent in reagents.items():
        db[tag] = with_local.eval(reagent)

    return db

def exec_app(src):
    return with_swmb.exec(src, get='app')

def eval_samples(src):
    return [
            make_namespace_recursive(x)
            for x in with_py.eval(src)
    ]

def eval_sample_group(src, default_key=None, default_member=None):
    default_key = default_key or {}
    default_member = default_member or {}

    schema = Schema({
        Optional('key', default={}): {str: with_swmb.eval},
        Optional('members', default=[]): [{str: with_swmb.eval}],
    })
    group = schema(src)
    group['key'] = {**default_key, **group['key']}

    key = make_namespace_recursive(group['key'])
    members = [
            make_namespace_recursive({**group['key'], **default_member, **x})
            for x in group['members'] or [group['key']]
    ]
    return stepwise_mol_bio.Group(key, members)

def match_protocol(protocol, expected, forbidden=[]):
    actual = protocol.format_text()

    prev = None
    print(actual.strip() + '\n')

    i = 0
    for x in always_iterable(expected):
        j = actual.find(x, i)

        if j == -1:
            print(f"Expected:\n  {x!r}\n")
            if prev:
                print(f"After:\n  {prev!r}\n")
            if isinstance(expected, str):
                delta = Differ().compare(
                        expected.splitlines(),
                        actual.splitlines(),
                )
                delta_str = '\n'.join(delta)
                print(f"Diff:\n{delta_str}")

            return False

        i = j + len(x)
        prev = x

    for x in always_iterable(forbidden):
        if x in actual:
            print(f"Didn't expect:\n  {x!r}")
            return False

    return True

def empty_ok(container):
    return Or(container, And('', lambda y: type(container)()))

def walk(f):
    def schema(xs):
        if isinstance(xs, list):
            return [f(x) for x in xs]
        if isinstance(xs, dict):
            return {k: f(v) for k, v in xs.items()}
        raise TypeError(f"expected list of dict, not {type(xs)}")
    return schema

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

    The factory function accepts arbitrary keyword arguments.  These arguments 
    are attached to the test function itself, and can be accessed from within 
    the test function via the request fixture: `request.node.function.kwargs`.
    """
    import inspect, sys
    from pathlib import Path

    def decorator(f):

        def factory(**meta_kwargs):
            frame = inspect.currentframe()
            try:
                ns = frame.f_back.f_globals

                test_func = copy_func(f)
                test_func.__module__ = Path(ns['__file__']).stem
                test_func.kwargs = meta_kwargs

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
        schema=[
            with_swmb.error_or('expected', 'forbidden'),
            defaults(forbidden=[]),
        ]
)
def test_protocol(app, expected, forbidden, error, request, disable_capture):
    app = exec_app(app)
    
    if not request.node.function.kwargs.get('disable_capture', False):
        disable_capture = nullcontext()

    with error:
        with disable_capture:
            protocol = app.protocol
        assert match_protocol(protocol, expected, forbidden)

@parametrize_from_file_factory(
        schema=defaults(stdout='^$', stderr='^$'),
)
@pytest.mark.slow
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout, stderr=stderr)

@parametrize_from_file_factory(
        schema=defaults(tags=[], forbidden=[]),
)
def test_freezerbox_make(db, tags, expected, forbidden, disable_capture):
    db = eval_db(db)
    tags = tags or list(db.keys())
    app = Make(db, tags)

    with disable_capture:
        protocol = app.protocol

    assert match_protocol(protocol, expected, forbidden)

@parametrize_from_file_factory(
        schema=defaults(errors={}),
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

