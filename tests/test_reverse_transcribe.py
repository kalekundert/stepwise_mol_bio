import stepwise
import stepwise_mol_bio.reverse_transcribe as rt
from param_helpers import *

@parametrize_from_file(
        schema=defaults(include_primers='True', kwargs={}),
)
def test_make_combos(group, reaction, include_primers, kwargs, expected):
    group = eval_sample_group(group)
    rxn = with_swmb.eval(reaction)
    kwargs = with_swmb.eval(kwargs)
    expected = with_py.eval(expected)

    assert rt.make_combos(group, rxn, **kwargs) == unordered(expected)

@parametrize_from_file
def test_setup_template(group, reaction, expected):
    group = eval_sample_group(group)
    rxn = with_swmb.eval(reaction)
    expected = with_swmb.eval(expected)

    rt.setup_template(group, rxn)
    assert rxn == expected

@parametrize_from_file
def test_setup_primer(group, reaction, expected):
    group = eval_sample_group(group)
    rxn = with_swmb.eval(reaction)
    expected = with_swmb.eval(expected)

    rt.setup_primer(group, rxn)
    assert rxn == expected

@parametrize_from_file(
        schema=with_swmb.error_or('expected'),
)
def test_denature_protocol(group, expected, error):
    group = eval_sample_group(group)
    with error:
        p = rt.plan_denature_protocol(group)
        assert match_protocol(p, expected)

@parametrize_from_file
def test_quick_reactions(group, expected):
    group = eval_sample_group(group)
    rxns = rt.plan_quick_reactions(group)
    assert match_protocol(rxns.protocol, expected)

@parametrize_from_file(
        schema=with_swmb.error_or('expected'),
)
def test_standard_protocol(group, expected, error):
    group = eval_sample_group(group)
    with error:
        p = rt.plan_standard_protocol(group)
        assert match_protocol(p, expected)

@parametrize_from_file
def test_incubation(group, expected):
    group = eval_sample_group(group)
    thermocycler = rt.plan_incubation(group)
    protocol = stepwise.Protocol() + thermocycler
    assert match_protocol(protocol, expected)

test_protocol()
test_cli()
