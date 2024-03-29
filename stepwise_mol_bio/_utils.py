#!/usr/bin/env python3

import sys
import byoc
import autoprop
import tidyexc
import freezerbox

from freezerbox import ReagentConfig, BaseProductConfig, iter_combo_makers
from byoc import Method, DocoptConfig
from appdirs import AppDirs
from inform import format_range, error
from more_itertools import all_equal, always_iterable, first, unique_everseen
from functools import partial
from pathlib import Path

app_dirs = AppDirs("stepwise_mol_bio")

@autoprop
class Main(byoc.App):
    usage_io = sys.stderr
    group_by = {}
    merge_by = {}

    @classmethod
    def main(cls):
        app = cls.from_bare()
        byoc.load(app, DocoptConfig)
        byoc.load(app, BaseProductConfig)
        
        try:
            app.protocol.print()
        except StepwiseMolBioError as err:
            error(err)

    @classmethod
    def make(cls, db, products, *, group_by=None, merge_by=None):
        if group_by is None:
            group_by = cls.group_by

        if merge_by is None:
            merge_by = cls.merge_by

        yield from iter_combo_makers(
                partial(cls._combo_maker_factory, db),
                map(cls._solo_maker_factory, products),
                group_by=group_by,
                merge_by=merge_by,
        )

    @classmethod
    def maker_from_reagent(cls, db, product):
        return cls._solo_maker_factory(product)

    @classmethod
    def protocols_from_makers(cls, makers):
        db = first(makers).db
        for maker in iter_combo_makers(
                    partial(cls._combo_maker_factory, db),
                    makers,
                    group_by=cls.group_by,
                    merge_by=cls.merge_by,
            ):
            yield maker.makers, maker.protocol


    def refresh(self):
        autoprop.clear_cache(self)

    def get_db(self):
        try:
            return self._db
        except AttributeError:
            self._db = freezerbox.load_db()
            return self._db

    def set_db(self, db):
        self._db = db

    @classmethod
    def _solo_maker_factory(cls, product):
        app = cls.from_bare()
        app.db = product.db
        app.products = [product]
        byoc.load(app, BaseProductConfig)
        return app

    @classmethod
    def _combo_maker_factory(cls, db):
        app = cls.from_bare()
        app.db = db
        # Shouldn't need to load `BaseProductConfig` here.  Any params that 
        # need it should be loaded in the solo makers and then copied into the 
        # combo maker.
        return app


class Cleanup(Main):

    product_tags = byoc.param(
            Method(lambda self: [x.tag for x in self.products], skip=AttributeError),
            default_factory=list,
    )

    def __bareinit__(self):
        super().__bareinit__()
        self.show_product_tags = False

    @classmethod
    def protocols_from_makers(cls, makers):
        db = first(makers).db
        combo_makers = list(iter_combo_makers(
                partial(cls._combo_maker_factory, db),
                makers,
                group_by=cls.group_by,
                merge_by=cls.merge_by,
        ))
        show_product_tags = (len(combo_makers) != 1)

        for maker in combo_makers:
            maker.show_product_tags = show_product_tags
            yield maker.makers, maker.protocol

@autoprop
class Bindable(metaclass=byoc.BareMeta):
    """
    Superclass for objects that can be bound to a Main/Cleanup instance, for 
    the purpose of gaining access to its FreezerBox database, BYOC config, 
    etc.

    See also: `bind()`
    """
    __config__ = []
    _use_app_configs = False

    @classmethod
    def partial(cls, *args, **kwargs):
        def factory(*args2, **kwargs2):
            return cls(*args, *args2, **{**kwargs, **kwargs2})
        return factory

    def __init__(self, *, db=None, **kwargs):
        self._db = db
        self._set_known_attrs(kwargs)

    def __init_subclass__(cls, **kwargs):
        cls._use_app_configs = kwargs.pop('use_app_configs', False)
        super().__init_subclass__(**kwargs)

    def bind(self, app, force=False):
        if not hasattr(self, 'app') or force:
            self.app = app
            self.on_bind(app, force=force)

    def on_bind(self, app, force=False):
        if self._use_app_configs:
            byoc.share_configs(app, self)

    def get_db(self):
        return self._db or self.app.db

    def set_db(self, db):
        self._db = db

    def _set_known_attrs(self, attrs):
        for attr, value in attrs.items():
            if attr not in self.__class__.__dict__:
                raise AttributeError(f"unknown attribute {attr!r}")
            setattr(self, attr, value)

@autoprop
class BindableReagent(Bindable):
    __config__ = [ReagentConfig]

    def __init__(self, tag, **kwargs):
        super().__init__(**kwargs)
        self.tag = tag

    def __str__(self):
        return self.tag

    def __repr__(self):
        return f'{self.__class__.__qualname__}({self.tag!r})'

    def __eq__(self, other):
        try:
            # Doesn't compare attributes, so be careful.
            return self.tag == other.tag
        except AttributeError:
            return NotImplemented

    def __hash__(self):
        return hash(self.tag)

@autoprop
class BindableReagents(Bindable):
    # Don't know if this is actually used anywhere, and should be deprecated 
    # anyways.

    def __init__(self, tags, **kwargs):
        super().__init__(**kwargs)
        self.tags = tags

    def __str__(self):
        return ','.join(self.tags)

    def __repr__(self):
        tag_reprs = ', '.join(map(repr, self.tags))
        return f'{self.__class__.__qualname__}({tag_reprs})'

    def __eq__(self, other):
        try:
            # Doesn't compare attributes, so be careful.
            return self.tags == other.tags
        except AttributeError:
            return NotImplemented



class StepwiseMolBioError(tidyexc.Error):
    pass

class ConfigError(StepwiseMolBioError):
    # For values that don't make sense, e.g. non-existent enzymes, etc.
    pass

class UsageError(StepwiseMolBioError):
    # For if the program isn't being used correctly, e.g. missing information.
    pass

def bind(app, bindables, iter=always_iterable, force=False):
    # Would `funcy.walk()` be a good way to handle any level of nesting?
    for bindable in iter(bindables):
        bindable.bind(app, force=force)
    return bindables

def try_except(expr, exc, failure, success=None):
    try:
        x = expr()
    except exc:
        return failure()
    else:
        return success() if success else x

def hanging_indent(text, prefix):
    from textwrap import indent
    if not isinstance(text, str):
        text = '\n'.join(map(str, text))
    if isinstance(prefix, int):
        prefix = ' ' * prefix
    return indent(text, prefix)[len(prefix):]

def merge_dicts(dicts):
    result = {}
    for dict in reversed(list(dicts)):
        result.update(dict)
    return result

def comma_list(x):
    return [x.strip() for x in x.split(',')]

def comma_tuple(x):
    return tuple(comma_list(x))

def comma_set(x):
    return set(comma_list(x))

def int_or_expr(x):
    return type_or_expr(int, x)

def float_or_expr(x):
    return type_or_expr(float, x)

def type_or_expr(type, x):
    if isinstance(x, str):
        return type(eval(x))
    else:
        return type(x)

def require_reagent(rxn, reagent):
    if reagent not in rxn:
        raise UsageError(f"reagent table missing {reagent!r}")

def merge_names(names):
    names = list(names)
    if all_equal(names):
        return names[0]

    name_str = ','.join(names)

    if len(name_str) > 30:
        name_str = ','.join(unique_everseen(names))
    if len(name_str) > 30:
        name_str = name_str[:29] + '…'

    return name_str

def match_len(x, n):
    # Something more generic than this might fit well in `more_itertools`.
    if isinstance(x, list):
        if len(x) != n:
            raise ValueError(f"expected {n} item(s), got {len(x)}")
        return x
    else:
        return n * [x]

def round_down_to_1_sig_fig(x):
    from math import log10, floor

    exponent = floor(log10(abs(x)))
    factor = 10**exponent
    mantissa = x / factor

    return floor(mantissa) * factor

def round_up_to_1_sig_fig(x):
    from math import log10, floor, ceil

    exponent = floor(log10(abs(x)))
    factor = 10**exponent
    mantissa = x / factor

    return ceil(mantissa) * factor

def format_sec(x):
    if x < 60:
        return f'{x}s'

    min = x // 60
    sec = x % 60

    return f'{min}m{f"{sec:02}" if sec else ""}'

def format_min(x):
    if x < 60:
        return f'{x}m'

    hr = x // 60
    min = x % 60

    return f'{hr}h{f"{sec:02}" if min else ""}'

