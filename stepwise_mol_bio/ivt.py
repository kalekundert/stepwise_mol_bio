#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import Key, DocoptConfig
from stepwise_mol_bio import Main, UsageError
from stepwise_mol_bio.gels.gel import Gel
from operator import not_
from inform import plural

@autoprop
class Ivt(Main):
    """\
Display a protocol for synthesizing in vitro transcribed (IVT) RNA using the
NEB HiScribe kit (E2040).

Usage:
    ivt <templates>... [options]

Arguments:
    <templates>
        The names to the DNA templates to transcribe.

Options:
    -d --dna-ng-uL CONC  [default: ${app.dna_ng_uL}]
        The concentration of the template DNA (in ng/µL).

    -D --dna-ng AMOUNT    [default: ${app.dna_ng}]
        How much template DNA to use (in ng).  NEB recommends 1 µg.  Lower 
        amounts will give proportionately lower yield, but will otherwise work 
        fine.  If the template is not concentrated enough to reach the given 
        quantity, the reaction will just contain as much DNA as possible 
        instead.

    -V --dna-uL VOLUME
        The volume of DNA to add to the reaction (in µL).  By default, this 
        will be inferred from --dna-ng and --dna-ng-uL.

    -s --short
        Indicate that all of the templates are shorter than 300 bp.  This 
        allows the reactions to be setup with less material, so that more 
        reactions can be performed with a single kit.

    -i --incubate HOURS
        How long to incubate the transcription reaction.  The default is 2h for 
        long transcripts and 4h for short transcripts (--short).

    -x --extra PERCENT      [default: ${app.extra_percent}]
        How much extra master mix to create.

    -R --no-rntp-mix
        Indicate that each you're not using a rNTP mix and that you need to add 
        each rNTP individually to the reaction.

    -c --cleanup METHOD     [default: ${app.cleanup}]
        Choose the method for removing free nucleotides from the RNA:

        none: Carry on the crude reaction mix.
        zymo: Zymo spin kits.
        ammonium: Ammonium acetate precipitation.

    -G --no-gel
        Don't include the gel electrophoresis step.
"""
    __config__ = [
            DocoptConfig(),
    ]

    templates = appcli.param(
            Key(DocoptConfig, '<templates>'),
    )
    _dna_uL = appcli.param(
            Key(DocoptConfig, '--dna-uL'),
            cast=float,
            default=None,
    )
    _dna_ng = appcli.param(
            Key(DocoptConfig, '--dna-ng'),
            cast=float,
            default=1000,
    )
    dna_ng_uL = appcli.param(
            Key(DocoptConfig, '--dna-ng-uL'),
            cast=float,
            default=500,
    )
    cleanup = appcli.param(
            Key(DocoptConfig, '--cleanup'),
            default='zymo',
    )
    short = appcli.param(
            Key(DocoptConfig, '--short'),
            default=False,
    )
    incubation_h = appcli.param(
            Key(DocoptConfig, '--incubate'),
            default=None,
    )
    rntp_mix = appcli.param(
            Key(DocoptConfig, '--no-rntp-mix', cast=not_),
            default=True,
    )
    gel = appcli.param(
            Key(DocoptConfig, '--no-gel', cast=not_),
            default=True,
    )
    extra_percent = appcli.param(
            Key(DocoptConfig, '--extra-percent'),
            cast=float,
            default=10,
    )

    def __init__(self, templates):
        self.templates = templates

    def get_reaction(self):
        ivt = stepwise.MasterMix("""\
                Reagent              Product #  Stock Conc      Volume  MM?
                ===================  =========  ==========  ==========  ===
                nuclease-free water                         to 20.0 μL  yes
                reaction buffer                        10x      2.0 μL  yes
                rATP                                100 mM      2.0 μL  yes
                rCTP                                100 mM      2.0 μL  yes
                rGTP                                100 mM      2.0 μL  yes
                rUTP                                100 mM      2.0 μL  yes
                HiScribe T7          NEB E2040         10x      2.0 μL  yes
                DNA template                     500 ng/μL      2.0 μL
        """)
        ivt.extra_percent = self.extra_percent
        ivt.num_reactions = len(self.templates)

        ivt['nuclease-free water'].order = 1
        ivt['reaction buffer'].order = 2
        ivt['rATP'].order = 3
        ivt['rCTP'].order = 3
        ivt['rGTP'].order = 3
        ivt['rUTP'].order = 3
        ivt['HiScribe T7'].order = 4
        ivt['DNA template'].order = 5

        if self.short:
            ivt['reaction buffer'].volume = '1.5 µL'
            ivt['rATP'].volume = '1.5 µL'
            ivt['rCTP'].volume = '1.5 µL'
            ivt['rGTP'].volume = '1.5 µL'
            ivt['rUTP'].volume = '1.5 µL'

        if self.rntp_mix:
            ivt['rNTP mix'].volume = 4 * ivt['rATP'].volume
            ivt['rNTP mix'].stock_conc = '100 mM'
            ivt['rNTP mix'].order = 3
            del ivt['rATP']
            del ivt['rCTP']
            del ivt['rGTP']
            del ivt['rUTP']

        ivt['DNA template'].name = ','.join(self.templates)
        ivt['DNA template'].stock_conc = self.dna_ng_uL, 'ng/µL'

        if self.dna_uL:
            ivt['DNA template'].volume = self.dna_uL, 'µL'
        else:
            ivt['DNA template'].volume = self.dna_ng / self.dna_ng_uL, 'µL'
        
        ivt.fix_volumes('DNA template', 'nuclease-free water')
        return ivt

    def get_protocol(self):
        p = stepwise.Protocol()

        ## Clean your bench
        p += """\
Wipe down your bench and anything you'll touch 
(pipets, racks, pens, etc.) with RNaseZap.
"""
        ## In vitro transcription
        ivt = self.reaction
        n = plural(ivt.num_reactions)
        p += f"""\
Setup {n:# in vitro transcription reaction/s} [1,2] by 
mixing the following reagents at room temperature 
in the order given:

{ivt}"""
        p.footnotes[1] = """\
https://tinyurl.com/y4a2j8w5
"""
        p.footnotes[2] = """\
I've found that T7 kits which have been in the 
freezer for more than ≈4 weeks seem to produce 
more degraded RNA.
"""
        t = self.incubation_h or (4 if self.short else 2)
        p += f"""\
Incubate at 37°C for {plural(t):# hour/s} [3].
"""
        p.footnotes[3] = """\
Use a thermocycler to prevent evaporation.
"""
        ## Purify product
        if self.cleanup == 'zymo':
            p += """\
Remove unincorporated ribonucleotides using Zymo 
RNA Clean & Concentrator 25 spin columns.
"""
        elif self.cleanup == 'ammonium':
            if self.short:
                raise UsageError("ammonium acetate precipitation cannot be used for short transcripts (<100 bp)")
            protocol += """\
Remove unincorporated ribonucleotides using
ammonium acetate precipitation:

- Add 1 volume (20 μL) 5M ammonium acetate to 
  each reaction.
- Incubate on ice for 15 min.
- Centrifuge at >10,000g for 15 min at 4°C.
- Wash pellet with 70% ethanol.
- Dissolve pellet in 20μL nuclease-free water.
"""
        elif self.cleanup == 'none':
            return p
        else:
            raise UsageError(f"unknown RNA clean-up method: '{self.cleanup}'")

#        ## Nanodrop concentration
#        p += """\
#Nanodrop to determine the RNA concentration.
#"""
#        ## Aliquot
#        p += """\
#Dilute (if desired) enough RNA to make several 
#10 μM aliquots and to run a gel.  Keep any left- 
#over RNA undiluted.  Flash-freeze in liquid N₂
#and store at -80°C.
#"""
#        ## Gel electrophoresis
#        if self.gel:
#            p += Gel('urea', len(self.templates)).protocol
#
        return p

    def get_dna_uL(self):
        return self._dna_uL

    def get_dna_ng(self):
        return self._dna_ng

    def set_dna_uL(self, dna_uL):
        self._dna_uL = dna_uL
        self._dna_ng = self.dna_ng_uL / dna_uL

    def set_dna_ng(self, dna_ng):
        self._dna_ng = dna_ng
        self._dna_uL = dna_ng / self.dna_ng_uL

    @appcli.on_load
    def _load_dna_quantity(self):
        if self._dna_uL:
            self.set_dna_uL(self._dna_uL)
        else:
            self.set_dna_ng(self._dna_ng)



if __name__ == '__main__':
    Ivt.main()

# vim: tw=50
