#!/usr/bin/env python3

"""\
Purify and concentrate nucleic acids by ethanol precipitation.

This protocol is primarily based on [Li2020].

Usage:
    ethanol_precipitation [<names>...] [options]

Arguments:
    <names>
        The names of the constructs to precipitate.

Options:
    -p --preset <name>                    [default: plasmid]
        There are four versions of the protocol, each optimized for a different 
        nucleic acid species.  Use this option to specify which version to use.  
        The names are case-insensitive:

        plasmid:
            Optimized with 10 kb circular plasmid.  This protocol is probably 
            also most appropriate for linear molecules of comparable size (e.g. 
            restriction digested plasmids).

        pcr:
            Optimized with 150 bp linear, doubled-stranded DNA.

        primer:
            Optimized with 20 nt single-stranded DNA.

        microrna:
            Optimized with 20 nt single-stranded RNA.

    -s --solvent <name>                     [default: etoh]
        The organic solvent to use for the precipitation.  The names are 
        case-insensitive.

        etoh: Ethanol
            - Gives higher yield than isopropanol for short RNA/DNA, and 
              comparable yield for longer DNA [Li2020].
            - Evaporates more easily after the precipitation.

        iproh: Isopropanol
            - Has been said to work better than ethanol for dilute samples, 
              although this was not tested by [Li2020].
            - Requires less volume, which may be beneficial when working with 
              large volumes.
            - Better at dissolving (and therefore removing) protein and 
              polysaccharide contaminants.
            - Precipitates more salt, resulting in higher salt contamination.

    -a --cation <name>
        The cation to use for the precipitation.  This is automatically 
        determined by the protocol, but you can specify a different choice 
        (e.g. based on what you have on hand).  The names are case-insensitive:

        na: Sodium acetate (NaAc)
        mg: Magnesium chloride (MgCl₂)

        Other cations were tested in [Li2020], but either NaAc or MgCl₂ was the 
        best in every condition.

    -c --carrier <name>
        The carrier, or coprecipitator, to add to the reaction.  This is 
        automatically determined by the protocol, but you can specify a 
        different choice (e.g. based on what you have on hand).  The names are 
        case-insensitive:

        lpa: linear polyacrylamide
            Not known to interfere with any downstream application.  Not 
            derived from a biological source, so very unlikely to have any 
            nucleic acid contamination.

        glycogen:
            Mostly inert, but may interfere with protein/DNA interactions 
            [Gaillard1990] and reverse transcription (at concentrations 
            >2 mg/mL).  Derived from biological source, so may contain trace 
            contaminating nucleic acids.  You can purchase glycogen crosslinked 
            to a blue dye, which makes the pellet even easier to see.

        tRNA: yeast tRNA
            Interferes with the quantification of the nucleic acid by Nanodrop, 
            which is problematic for many applications.

    -b --buffer <name>                      [default: water]
        The aqueous buffer to resuspend the preciptated nucleic acid in.

    -v --buffer-volume <µL>
        The volume of resuspension buffer to use, in µL

    -I --no-incubation
        Exclude the incubation step.

    -W --no-wash
        Exclude the wash step.
            
References:
    Li Y et al.  A systematic investigation of key factors of nucleic acid 
    precipitation toward optimized DNA/RNA isolation.  BioTechniques 68, 
    191–199 (2020).

    Gaillard C, Strauss F.  Ethanol precipitation of DNA with linear 
    polyacrylamide as carrier.  Nucleic Acids Res. 18(2), 378 (1990).

    Sambrook J & Russell DW.  Standard ethanol precipitation of DNA in 
    microcentrifuge tubes.  Cold Spring Harb Protoc (2006).
"""

import stepwise
import autoprop
import textwrap
from inform import plural
from fractions import Fraction
from stepwise_mol_bio import Main

@autoprop
class EthanolPrecipitation(Main):
    presets = {
            'plasmid': {
                'solvent': 'etoh',
                'solvent_volume': {
                    'etoh': 3,
                    'iproh': 1,
                },

                'cation': {
                    'etoh': 'na',
                    'iproh': 'na',
                },

                'carrier': {
                    'etoh': 'lpa',
                    'iproh': 'lpa',
                },

                'incubation_time': None,
                'incubation_temp_C': None,

                'centrifugation_time_min': 60,
                'centrifugation_temp_C': 4,
                'centrifugation_speed': '>7500g',
            },
            'pcr': {
                'solvent': 'etoh',
                'solvent_volume': {
                    'etoh': 2,
                    'iproh': Fraction(3,4),
                },

                'cation': {
                    'etoh': 'mg',
                    'iproh': 'mg',
                },

                'carrier': {
                    'etoh': 'glycogen',
                    'iproh': 'lpa',
                },

                'incubation_time': 'overnight',
                'incubation_temp_C': -20,

                'centrifugation_time_min': 60,
                'centrifugation_temp_C': 4,
                'centrifugation_speed': '>7500g'
            },
            'primer': {
                'solvent': 'etoh',
                'solvent_volume': {
                    'etoh': 4,
                    'iproh': 1,
                },

                'cation': {
                    'etoh': 'na',
                    'iproh': 'mg',
                },

                'carrier': {
                    'etoh': 'glycogen',
                    'iproh': 'glycogen',
                },

                'incubation_time': 'overnight',
                'incubation_temp_C': 4,

                'centrifugation_time_min': 60,
                'centrifugation_temp_C': 4,
                'centrifugation_speed': '>18000g',
            },
            'microrna': {
                'solvent': 'etoh',
                'solvent_volume': {
                    'etoh': 4,
                    'iproh': Fraction(3,4),
                },

                'cation': {
                    'etoh': 'mg',
                    'iproh': 'na',
                },

                'carrier': {
                    'etoh': 'glycogen',
                    'iproh': 'lpa',
                },

                'incubation_time': 'overnight',
                'incubation_temp_C': -20,

                'centrifugation_time_min': 60,
                'centrifugation_temp_C': 4,
                'centrifugation_speed': '>21000',
            },
    }
    solvents = {
            'etoh': {
                'name': '100% ethanol',
            },
            'iproh': {
                'name': 'isopropanol',
            },
    }
    carriers = {
            'trna': {
                'name': "yeast tRNA",
                'conc': "20 ng/µL",
            },
            'glycogen': {
                'name': "glycogen",
                'conc': "50 ng/µL",
            },
            'lpa': {
                'name': "linear polyacrylamide (LPA)",
                'conc': "20 ng/µL",
            },
    }
    cations = {
            'na': {
                'name': "sodium acetate, pH=5.2",
                'conc': "300 mM",
            },
            'mg': {
                'name': "magnesium chloride (MgCl₂)",
                'conc': "10 mM",
            },
    }

    def __init__(self, preset):
        self.preset = preset
        self.names = None
        self._solvent = None
        self._cation = None
        self._carrier = None
        self.buffer = None
        self.buffer_volume_uL = None
        self.incubation = True
        self.wash = True

    @classmethod
    def from_docopt(cls, args):
        self = cls(args['--preset'])
        self.names = args['<names>']
        self.solvent = args['--solvent']
        self.cation = args['--cation']
        self.carrier = args['--carrier']
        self.wash = not args['--no-wash']
        self.buffer = args['--buffer']
        self.buffer_volume_uL = args['--buffer-volume']
        return self

    def get_protocol(self):
        d = self.defaults
        p = stepwise.Protocol()

        if self.names:
            s = stepwise.Step(f"Purify {','.join(names)} by ethanol precipitation [1]:")
        else:
            s = stepwise.Step("Perform an ethanol precipitation [1]:")

        p += s

        s += f"""\
                Add {self.cation_name} to {self.cation_conc}."""
        s += f"""\
                Add {self.carrier_name} to {self.carrier_conc}."""
        s += f"""\
                Add {plural(self.solvent_volume):# volume/s} 
                {self.solvent_name} and mix well."""
        s += f"""\
                If necessary, divide the sample between microfuge tubes
                such that none holds more than 400 µL."""

        if self.incubation and (t := self.incubation_time):
            incubation_time = "overnight" if t == 'overnight' else f"for {t}"
            s += f"""\
                    Incubate at {self.incubation_temp_C}°C {incubation_time} 
                    [2]."""

        s += f"""\
                Centrifuge {self.centrifugation_speed}, 
                {self.centrifugation_time_min} min, 
                {self.centrifugation_temp_C}°C [3].  Remove the supernatant, 
                but save it in case the precipitation needs to be repeated."""

        if self.wash:
            s += f"""\
                    Add 800 µL 70% ethanol."""
            s += f"""\
                    Centrifuge {self.centrifugation_speed}, 2 min, 
                    {self.centrifugation_temp_C}°C and discard supernatant."""
        s += """\
                Leave the tube open at room temperature until ethanol has 
                evaporated [4]."""

        s += f"""\
                Resuspend the pellet in {f'{self.buffer_volume_uL} µL' or
                'any volume'} of {self.buffer} [5]."""

        p.footnotes[1] = textwrap.dedent("""\
                Li2020: 10.2144/btn-2019-0109
                Sambrook2006: 10.1101/pdb.prot4456
        """)
        p.footnotes[2] = stepwise.Footnote("""\
                DNA can be stored indefinitely in ethanolic solutions at either 
                0°C or −20°C.
        """)
        p.footnotes[3] = stepwise.Footnote("""\
                This protocol was optimized for 100 ng/µL nucleic acid.  If 
                your sample is substantially more dilute, it may be necessary 
                to compensate by increasing the centrifugation time and/or 
                speed.
        """)
        p.footnotes[4] = stepwise.Footnote("""\
                Do not dry pellets of nucleic acid in a lyophilizer, as this 
                causes denaturation of small (<400-nucleotide) fragments of DNA 
                and greatly reduces the recovery of larger fragments of DNA. 

                If necessary, the open tube containing the redissolved DNA can 
                be incubated for 2-3 minutes at 45°C in a heating block to 
                allow any traces of ethanol to evaporate.
        """)
        p.footnotes[5] = stepwise.Footnote("""\
                Up to 50% of the DNA is smeared on the wall of the tube. To 
                recover all of the DNA, push a bead of fluid backward and 
                forward over the appropriate quadrant of wall with a pipette 
                tip.
        """)

        p.prune_footnotes()
        return p

    def get_defaults(self):
        return self.presets[self.preset]

    def get_solvent(self):
        return (self._solvent or self.defaults['solvent']).lower()

    def set_solvent(self, solvent):
        self._solvent = solvent

    def get_solvent_name(self):
        return self.solvents[self.solvent]['name']

    def get_solvent_volume(self):
        return self.defaults['solvent_volume'][self.solvent]

    def get_cation(self):
        return (self._cation or self.defaults['cation'][self.solvent]).lower()

    def set_cation(self, cation):
        self._cation = cation

    def get_cation_name(self):
        return self.cations[self.cation]['name']

    def get_cation_conc(self):
        return self.cations[self.cation]['conc']

    def get_carrier(self):
        return (self._carrier or self.defaults['carrier'][self.solvent]).lower()

    def set_carrier(self, carrier):
        self._carrier = carrier

    def get_carrier_name(self):
        return self.carriers[self.carrier]['name']

    def get_carrier_conc(self):
        return self.carriers[self.carrier]['conc']

    def get_incubation_temp_C(self):
        return self.defaults['incubation_temp_C']

    def get_incubation_time(self):
        return self.defaults['incubation_time']

    def get_centrifugation_speed(self):
        return self.defaults['centrifugation_speed']

    def get_centrifugation_time_min(self):
        return self.defaults['centrifugation_time_min']

    def get_centrifugation_temp_C(self):
        return self.defaults['centrifugation_temp_C']

if __name__ == '__main__':
    EthanolPrecipitation.main(__doc__)
