#!/usr/bin/env python3

import stepwise
import stepwise_mol_bio
import freezerbox
import parametrize_from_file

from voluptuous import Schema, Invalid, Coerce, And, Or, Optional
from parametrize_from_file.voluptuous import Namespace, empty_ok
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
    for x in expected:
        j = actual.find(x, i)

        if j == -1:
            print(f"Expected:\n  {x!r}")
            if prev:
                print(f"After:\n  {prev!r}")

            return False

        i = j + len(x)
        prev = x

    for x in forbidden:
        if x in actual:
            print(f"Didn't expect:\n  {x!r}")
            return False

    return True

def noop(x):
    return x

