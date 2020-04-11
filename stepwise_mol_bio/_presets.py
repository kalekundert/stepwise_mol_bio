#!/usr/bin/env python3

import stepwise

class Presets:

    def __init__(self, presets):
        self.presets = presets

    @classmethod
    def from_config(cls, dotted_key):
        try:
            config = stepwise.load_config()
            config = _get_dotted_key(config, dotted_key)
            return cls(config.data)

        except KeyError:
            return cls({})

    def load(self, key):
        """
        Return a dictionary containing all of the parameters---including 
        inherited parameters---associated with the specified preset.
        """
        return _load_preset(self.presets, key)

    def format_briefs(self, default=''):
        rows = []
        for key in self.presets:
            preset = self.load(key)
            rows.append([
                    f'{key}:',
                    preset.get('brief', default.format(**preset)),
            ])
        return stepwise.tabulate(rows, truncate='-x', max_width=71)

def _load_preset(presets, key):
    from inform import did_you_mean
    from ._utils import ConfigError

    try:
        preset = presets[key]
    except KeyError:
        raise ConfigError(f"no preset {key!r}, did you mean {did_you_mean(str(key), presets)!r}") from None

    if 'inherit' not in preset:
        return preset
    else:
        parent = _load_preset(presets, preset['inherit'])
        merged = {**parent, **preset}; del merged['inherit']
        return merged

def _get_dotted_key(dict, dotted_key):
        for key in dotted_key.split('.'):
            dict = dict[key]
        return dict


