#!/usr/bin/env python3

"""\
Perform restriction digests using the protocol recommended by NEB.

Usage:
    digest <templates> <enzyme> [-d <ng>] [-D <ng/µL>] [-v <µL>] [-g]

Arguments:
    <templates>
        The DNA to digest.  Use commas to specify multiple templates.

    <enzyme>
        The restriction enzyme to use.  Only NEB enzymes are currently 
        supported.  If you are using an "HF" enzyme, specify that explicitly.  
        For example, "HindIII" and "HindIII-HF" have different protocols.  
        Enzyme names are case-insensitive.

Options:
    -d --dna <µg>               [default: 1]
        The amount of DNA to digest, in µg.

    -D --dna-stock <ng/µL>      [default: 200]
        The stock concentration of the DNA template, in ng/µL.

    -v --target-volume <µL>     [default: 10]
        The ideal volume for the digestion reaction.  Note that the actual 
        reaction volume may be increased to ensure that the volume of enzyme 
        (which is determined by the amount of DNA to digest, see --dna) is less 
        than 10% of the total reaction volume, as recommended by NEB.

    -g --genomic
        Indicate that genomic DNA is being digested.  This will double the 
        amount of enzyme used, as recommended by NEB.
"""

import docopt
import json
import requests
import stepwise
import autoprop
from pathlib import Path
from inform import Error, plural, did_you_mean
from stepwise_mol_bio.utils import app

@autoprop
class RestrictionDigest:

    def __init__(self, templates, enzyme):
        self.templates = templates
        self.enzyme = enzyme

        self.dna_ug = 1
        self.dna_stock_nguL = 200
        self.target_volume_uL = 10
        self.is_genomic = False

    def get_reaction(self):
        # Define a prototypical restriction digest reaction.  Stock 
        # concentrations for BSA, SAM, and ATP come from the given catalog 
        # numbers.

        rxn = stepwise.MasterMix.from_text("""\
                Reagent  Catalog      Stock    Volume  MM?
                =======  =======  =========  ========  ===
                water                        to 50 µL  yes
                DNA               200 ng/µL      5 µL
                BSA        B9000   20 mg/mL   0.25 µL  yes
                SAM        B9003      32 mM   0.25 µL  yes
                ATP        P0756      10 mM   0.25 µL  yes
                buffer                  10x      5 µL  yes
                enzyme              10 U/µL      1 µL  yes
        """)

        # Plug in the parameters the user requested.

        rxn.num_reactions = len(self.templates)

        rxn['DNA'].name = ','.join(self.templates)
        rxn['DNA'].hold_conc.stock_conc = self.dna_stock_nguL, 'ng/µL'
        
        rxn['enzyme'].name = self.enzyme['name']
        rxn['enzyme'].hold_conc.stock_conc = self.enzyme['concentration'] / 1000, 'U/µL'

        if self.is_genomic:
            rxn['enzyme'].volume *= 2

        rxn['buffer'].name = self.enzyme['recommBuffer']
        
        if not self.enzyme['supplement']['bsa']:
            del rxn['BSA']

        if conc_uM := self.enzyme['supplement']['sam']:
            rxn['SAM'].hold_stock.conc = conc_uM / 1000, 'mM'
        else:
            del rxn['SAM']

        if conc_mM := self.enzyme['supplement']['atp']:
            rxn['ATP'].hold_stock.conc = conc_mM, 'mM'
        else:
            del rxn['ATP']

        # Update the reaction volume.  This takes some care, because the 
        # reaction volume depends on the enzyme volume, which in turn depends 
        # on the DNA quantity.

        k = self.dna_ug / 1  # The prototype reaction has 1 µg DNA.
        dna_vol = k * rxn['DNA'].volume
        enz_vol = k * rxn['enzyme'].volume

        rxn.hold_ratios.volume = max(
                stepwise.Quantity(self.target_volume_uL, 'µL'),
                10 * enz_vol,

                # This is a bit of a hack.  The goal is to keep the water 
                # volume non-negative, but it won't necessarily work if there 
                # are supplements.
                10/9 * (dna_vol + enz_vol),
        )
        rxn['DNA'].volume = dna_vol
        rxn['enzyme'].volume = enz_vol

        return rxn

    def get_protocol(self):
        protocol = stepwise.Protocol()
        rxn = self.reaction

        protocol += f"""\
Setup {plural(rxn.num_reactions):# {rxn['enzyme'].name} digestion/s} [1]:

{rxn}
"""
        protocol += f"""\
Incubate at the following temperatures [2]:

- {self.enzyme['incubateTemp']}°C for {'5–15 min' if self.enzyme['timeSaver'] else '1 hour'}.
- {self.enzyme['heatInactivationTemp']}°C for {self.enzyme['heatInactivationTime']} min.
"""
        protocol.footnotes[1] = """\
NEB recommends 5–10 units of enzyme per µg DNA 
(10–20 units for genomic DNA).  Enzyme volume 
should not exceed 10% of the total reaction 
volume to prevent star activity due to excess 
glycerol.
"""
        protocol.footnotes[2] = """\
The heat inactivation step is not necessary if 
the DNA will be purified before use.
"""
        return protocol


class RestrictionDigestError(Error):
    pass

class CantDownloadEnzymes(RestrictionDigestError):
    template = """\
            Failed to download restriction enzyme data from NEB.  Make sure the 
            internet is connected and that the following URL is reachable:  

            {url}

            Note that the data only needs to be downloaded once.  After that, 
            this script should work offline."""
    wrap = True

    def __init__(self, url):
        super().__init__(url=url)


class UnknownEnzyme(RestrictionDigestError):
    template = "No such enzyme {enzyme_name!r}.  Did you mean {did_you_mean!r}?"
    wrap = True

    def __init__(self, enzyme_name, enzymes):
        super().__init__(
                enzyme_name=enzyme_name,
                enzymes=enzymes,
                did_you_mean=did_you_mean(enzyme_name, enzymes),
        )

def load_neb_enzyme(name):
    enzymes = load_neb_enzymes()
    enzymes_lower = {k.lower(): v for k, v in enzymes.items()}

    try:
        return enzymes_lower[name.lower()]
    except KeyError:
        raise UnknownEnzyme(name, enzymes)

def load_neb_enzymes():
    cache = Path(app.user_cache_dir) / 'neb' / 'restriction_enzymes.json'
    cache.parent.mkdir(parents=True, exist_ok=True)

    try:
        url = 'http://nebcloner.neb.com/data/reprop.json'
        data = requests.get(url).json()

        with cache.open('w') as f:
            json.dump(data, f)

    except requests.exceptions.ConnectionError:
        if not cache.exists():
            raise NoEnzymeData(url)

        with cache.open() as f:
            data = json.load(f)

    return data


if __name__ == '__main__':
    try:
        args = docopt.docopt(__doc__)
        templates = [x.strip() for x in args['<templates>'].split(',')]
        enzyme = load_neb_enzyme(args['<enzyme>'])

        digest = RestrictionDigest(templates, enzyme)
        digest.dna_ug = float(args['--dna'])
        digest.dna_stock_nguL = float(args['--dna-stock'])
        digest.target_volume_uL = float(args['--target-volume'])
        digest.is_genomic = args['--genomic']

        print(digest.protocol)

    except Error as err:
        err.report()
