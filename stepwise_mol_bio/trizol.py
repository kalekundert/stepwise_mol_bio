#!/usr/bin/env python3

import stepwise
import autoprop
import byoc

from stepwise_mol_bio import Main, UsageError
from stepwise import pl, ul
from byoc import DocoptConfig, Key
from functools import partial

@autoprop
class Trizol(Main):
    """
Extract total RNA using TRIzol, i.e. guanidinium thiocyanate-phenol-chloroform extraction.

Usage:
    trizol [-v <mL>] [-e <µL>] [-s <type>] [--high-fat-content]

Options:
    -v --trizol-volume <mL>
        How much TRIzol reagent to use, in mL.  If not specified, the protocol 
        will recommend the amount of reagent to use based on the amount of 
        sample you have.  If specified, the protocol will give exact volumes 
        for each step, saving yourself some mental math.

    -e --elution-volume <µL>
        What volume, in µL, to resuspend the RNA in after the isopropanol 
        precipitation.

    -s --sample-type <name>
        What kind of sample you are extracting RNA from.  This can be either 
        "tissue", "monolayer", or "suspension".  This only affects how the cell 
        lysis is performed.

    --high-fat-content
        Indicate if the sample has high fat content.  If so, an extra centrifugation step will be included to clarify the sample.
"""
    __config__ = [DocoptConfig]

    def __init__(self, trizol_volume_mL=None):
        if trizol_volume_mL:
            self.trizol_volume_mL = trizol_volume_mL

    def get_protocol(self):
        p = stepwise.Protocol()
        p += self.lysis_protocol
        p += self.extraction_protocol
        p += self.precipitation_protocol
        return p

    def get_lysis_protocol(self):
        p = stepwise.Protocol()
        p += pl(
                f"Lyse and homogenize samples in TRIzol reagent{p.add_footnotes(self.url)}:",
                steps := ul(),
        )

        v = self.trizol_volume_mL

        if self.sample_type == 'tissue':
            if v: steps += f"Add {v:g} mL TRIzol reagent."
            else: steps += "Add 1 mL TRIzol reagent per 50-100 mg tissue."
            steps += "Homogenize using a homogenizer."

        elif self.sample_type == 'monolayer':
            steps += "Remove growth media."
            if v: steps += f"Add {v:g} mL TRIzol reagent directly to the culture dish."
            else: steps += "Add 300-400 µL TRIzol reagent per 10⁵-10⁷ cells directly to the culture dish."
            steps += "Pipet the lysate up and down several times to homogenize."

        elif self.sample_type == 'suspension':
            steps += "Collect cells by centrifugation, discard supernatant."
            if v: steps += f"Add {v:g} mL TRIzol reagent."
            else: steps += "Add 750 µL TRIzol reagent per 250 µL sample."
            steps += "Pipet the lysate up and down several times to homogenize."

        else:
            raise UsageError(f"unknown sample type: {self.sample_type!r}")

        if self.high_fat_content:
            steps += "Centrifuge the lysate for 5 min at 12000g and 4-10°C."
            steps += "Transfer the clear supernatant to a clean tube."

        steps += "Incubate for 5 min at room temperature."

        return p

    def get_extraction_protocol(self):
        uL_per_mL = partial(
                volume_per_mL_trizol,
                trizol_mL=self.trizol_volume_mL,
                volume_unit='µL',
        )

        p = stepwise.Protocol()
        p += pl(
                f"Extract RNA from TRIzol reagent{p.add_footnotes(self.url)}:",
                ul(
                    f"Add {uL_per_mL('chloroform', 200)}.",
                    "Vortex vigorously.",
                    "Incubate for 2-3 min at room temperature.",
                    "Centrifuge for 15 min at 12000g and 4°C.",
                    f"Transfer the aqueous phase (top, not pink, ≈{uL_per_mL(None, 500)}) for each sample to a clean tube, taking care to avoid transferring any of the organic phase.",
                ),
        )
        return p

    def get_precipitation_protocol(self):
        uL_per_mL = partial(
                volume_per_mL_trizol,
                trizol_mL=self.trizol_volume_mL,
                volume_unit='µL',
        )
        mL_per_mL = partial(
                volume_per_mL_trizol,
                trizol_mL=self.trizol_volume_mL,
                volume_unit='mL',
        )

        p = stepwise.Protocol()
        p += pl(
                f"Concentrate and purify the RNA by isopropanol precipitation{p.add_footnotes(self.url)}:",
                ul(
                    "Add 5-10 µg RNAse-free glycogen.",
                    f"Add {uL_per_mL('isopropanol', 500)}.",
                    "Incubate at 4°C for 10 min."
                ),
                ul(
                    "Centrifuge for 10 min at 12,000g and 4°C.",
                    "Carefully pipet off all supernatant.",
                ),
                ul(
                    f"Resuspend pellet in {mL_per_mL('75% ethanol', 1)}.",
                    "Vortex briefly.",
                    "Centrifuge for 5 min at 7,500g and 4°C.",
                ),
                ul(
                    "Carefully pipet off all supernatant.",
                    "Air dry for 5-10 min.",
                    f"Resuspend RNA in {self.elution_volume_uL} µL water.",
                    "Incubate at 55-60°C for 10-15 min.",
                ),
        )
        return p

    trizol_volume_mL = byoc.param(
            Key(DocoptConfig, '--trizol-volume', cast=float),
            default=None,
    )
    elution_volume_uL = byoc.param(
            Key(DocoptConfig, '--elution-volume'),
            default='20-50',
    )
    sample_type = byoc.param(
            Key(DocoptConfig, '--sample-type'),
            default='suspension',
    )
    high_fat_content = byoc.param(
            Key(DocoptConfig, '--high-fat-content'),
            default=False,
    )
    url = 'https://tinyurl.com/yc6es8av'

def volume_per_mL_trizol(reagent, reagent_volume, volume_unit, trizol_mL):
    if trizol_mL:
        parts = [f'{reagent_volume * trizol_mL:.0f} {volume_unit}', reagent]
    else:
        parts = [f'{reagent_volume:.0f} {volume_unit}', reagent, 'per 1 mL TRIzol reagent']

    return ' '.join(x for x in parts if x)


if __name__ == '__main__':
    Trizol.entry_point()
