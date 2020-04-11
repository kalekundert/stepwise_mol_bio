#!/usr/bin/env python3

from stepwise_mol_bio import _presets as module
from utils import parametrize_via_toml

@parametrize_via_toml('test_presets.toml')
def test_load_preset(presets, key, expected):
    assert module._load_preset(presets, key) == expected

@parametrize_via_toml('test_presets.toml')
def test_get_dotted_key(dict, dotted_key, expected):
    assert module._get_dotted_key(dict, dotted_key) == expected


