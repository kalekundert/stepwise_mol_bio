#!/usr/bin/env python3

import byoc
import sys, inspect, functools

from . import StepwiseMolBioError
from byoc import DocoptConfig
from decopatch import function_decorator, DECORATED
from inform import error
from more_itertools import one, first
from types import SimpleNamespace
from freezerbox import group_by_identity

# Group object:
# - easy to pass around
# - clear about common vs variable
# - can't use generic tools

# Sample class:
# - Constructor that can initialize known attrs.
#   - Bindable already does this, but might want to move that feature.
# - Shortcut to make factory with partial args.
#   - Basically: take kwargs, invoke ctor with those args.

class App(byoc.App):
    usage_io = sys.stderr

    def main(self):
        byoc.load(self, DocoptConfig)

        self.protocol.print()
        # try:
        #     self.protocol.print()
        # except StepwiseMolBioError as err:
        #     error(err)

    @property
    def usage_vars(self):
        # Most of the BYOC configs are associated with the Sample class, but we 
        # would usually like to incorporate some information from those configs 
        # into the usage text (e.g. paths to the config files).  We accommodate 
        # this by instantiating a sample object for just this purpose.
        return {'sample': self.Sample.from_bare()}

class Group:

    def __init__(self, key, members, parent=None):
        self.__key = key
        self.__parent = parent
        self.__members = members

    def __iter__(self):
        yield from self.__members

    def __getattr__(self, attr):
        try:
            return getattr(self.__key, attr)
        except AttributeError:
            if self.__parent is not None:
                return getattr(self.__parent, attr)
            else:
                raise

@function_decorator
def group_samples(*attrs, f=DECORATED):
    """
    Allow a function to be called with either samples or groups, whichever is 
    convenient for the caller.

    The function itself will always be passed a group object, since that makes 
    it easy to access shared attributes. The function must either (i) only take 
    one argument or (ii) have exactly one argument annotated as a `Group`.  The 
    value associated with that argument will be extracted, converted to a group 
    if necessary, then passed to the function.

    It is not necessary to specify any *attrs*, but if you do, the decorated 
    function will also have a `for_groups_in` generator attached to it.  This 
    generator will have the same signature as the original function, except 
    that the group argument will now be expected to be an iterable.  That 
    iterable will be broken into groups according to the given attributes, then 
    each group will be passed to the original function and yielded.
    """

    class NotFound(Exception):
        pass

    def find_group(f, args, kwargs):
        sig = inspect.signature(f)
        try:
            group_param = one(
                    (k for k, v in sig.parameters.items()
                    if v.annotation is Group),
                    too_short=NotFound,
                    too_long=TypeError("at most one parameter can be annotated as Group"),
            )
        except NotFound:
            group_param = first(sig.parameters)

        bound_args = sig.bind(*args, **kwargs)
        group_arg = bound_args.arguments[group_param]

        def call_with_group(group):
            bound_args.arguments[group_param] = group
            return f(*bound_args.args, **bound_args.kwargs)

        return group_arg, call_with_group

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        group, call_with_group = find_group(f, args, kwargs)
        if not isinstance(group, Group):
            group = Group(group, [group])
        return call_with_group(group)

    def for_group_in(*args, **kwargs):
        group, call_with_group = find_group(f, args, kwargs)
        for subgroup in group_by_attrs(group, f.attrs):
            yield call_with_group(subgroup)

    # Automatically determine attrs
    # =============================
    # It would be nice to determine the *attrs* argument automatically, since 
    # in principle all the information is there in the code.  There are a few 
    # ways to try doing this:
    #
    # - Call the function with a mock object and record attribute accesses.
    #
    #   This probably isn't a very good idea.  There might be conditionals that 
    #   the code doesn't visit, and errors could be triggered if there are 
    #   checks for certain attribute values.  That said, it would just be a 
    #   default.
    #
    # - Parse the AST:
    #
    #   Getting an AST requires having access to the source code, which 
    #   function objects don't really have (they can just try to look it up in 
    #   the filesystem knowing their qualified names).  This means the AST is 
    #   unavailable for C-extensions and code defined in the interpreter.  That 
    #   said, source should be available for any application I can envision.
    # 
    #   Pseudocode:
    #   - Look for all uses of the group argument
    #     - If attribute access: record that
    #     - If assignment: start looking for that variable too
    #     - If function call: recurse into that function
    #
    # - Parse the bytecode:
    #
    #   Bytecode is lower-level than the AST, but it's available for any 
    #   function object.  There are two things I'd need to look for: attribute 
    #   accesses and function calls.  Attribute access aren't too bad.  It 
    #   seems like they'd always be LOAD_FAST(0) followed by LOAD_ATTR.  (The 0 
    #   is because the group is always the first argument.)  Function calls are 
    #   more difficult, though.  The pattern will be LOAD_GLOBAL (for the 
    #   function), followed by LOAD_FAST for each argument, followed by 
    #   CALL_FUNCTION.  There are also the CALL_FUNCTION_KW and 
    #   CALL_FUNCTION_EX opcodes, with slightly different usages.  It probably 
    #   could be done, but it'd be intricate.

    wrapper.for_group_in = for_group_in
    wrapper.attrs = f.attrs = set(attrs)

    return wrapper

def group_by_attrs(items, attrs):

    def by_attrs(item):
        return {
                k: getattr(item, k)
                for k in attrs
        }

    parent = item if isinstance(items, Group) else None

    for key, group in group_by_identity(items, key=by_attrs):
        yield Group(SimpleNamespace(**key), group, parent=parent)

