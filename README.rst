****************************
Stepwise — Molecular Biology
****************************

.. image:: https://img.shields.io/pypi/v/stepwise_mol_bio.svg
   :alt: Last release
   :target: https://pypi.python.org/pypi/stepwise_mol_bio

.. image:: https://img.shields.io/pypi/pyversions/stepwise_mol_bio.svg
   :alt: Python version
   :target: https://pypi.python.org/pypi/stepwise_mol_bio

..
  .. image:: https://img.shields.io/readthedocs/stepwise_mol_bio.svg
     :alt: Documentation
     :target: https://stepwise_mol_bio.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/github/actions/workflow/status/kalekundert/stepwise_mol_bio/test_and_release.yml?branch=master
   :alt: Test status
   :target: https://github.com/kalekundert/stepwise_mol_bio/actions

.. image:: https://img.shields.io/coveralls/kalekundert/stepwise_mol_bio.svg
   :alt: Test coverage
   :target: https://coveralls.io/github/kalekundert/stepwise_mol_bio?branch=master

.. image:: https://img.shields.io/github/last-commit/kalekundert/stepwise_mol_bio?logo=github
   :alt: Last commit
   :target: https://github.com/kalekundert/stepwise_mol_bio

Installation
============
Install ``stepwise_mol_bio`` using ``pip``::

    $ pip install stepwise_mol_bio

Writing a protocol
==================
While there are a lot of tools to help make writing protocols easier, it's not 
necessarily easy to see how they all fit together.  Below are a set of 
guidelines to help get you on the right track:

- The protocols in this repository are much more sophisticated than the usual 
  protocols you'd write.  That's because these protocols represent very common 
  techniques, and are meant to be *very* reuseable.

- Every protocol should:

  - Read options from the command-line.  Many protocols support the following 
    options:

    - -p, --preset
    - -v, --volume
    - -c, --conc
    - -n, --num-reactions
    - -u, --product

  - Read options from any stepwise configuration file.
  - Read options from the FreezerBox database, if available.
  - Be completely configurable from python

- The typical architecture for a script:

  - Use `appcli` to read attributes from the various data sources mentioned 
    above.

  - Typical classes:

    - `App`:

      - Inherit from `Main` or `Cleanup`.
      - Represents the complete protocol, which is usually a single master mix 
        and some follow-up steps.
      - The command-line and FreezerBox interfaces will be derived from this 
        object.

    - `Reaction`:

      - Inherit from `Argument`.
      - Represents a single reaction.  The app will hold a list of reactions.
      - There's a bit of an art to deciding which attributes go in the app vs 
        the reaction.  In general, anything that could vary between reactions 
        (e.g. stock concentrations) should go in the reaction.  Anything that 
        can't vary without it becoming impossible to make a master mix (e.g.  
        final concentrations, volumes) should go in the app.

    - `Reagent`:

      - Inherit from `Argument`.
      - Represents the individual reagents that make up the reaction, e.g.  
        template and primers for PCR, DNA and enzyme for restriction digests, 
        backbone and fragments for assemblies, etc.

    - The reaction and reagent classes sometimes have different names, if 
      there's a term that makes more sense for the specific protocol at hand.

      - For example, the PCR protocol uses the term "amplicon" instead of 
        "reaction", and has reagents named "template" and "primer".

    - Sometimes it makes sense to combine the reaction and reagent classes:

      - For example, the digest protocol would naturally have a reaction class 
        with enzyme and template attributes.  But the enzymes are intrinsic to 
        the master mix (e.g. two digestions with different enzymes would just 
        be done by invoking the digest command twice), and therefore are stored 
        in the app itself.  That leaves only the template, which on its own is 
        more like a reagent.

  - The reactions/reagents must be bound to the app.  This is how the reagents 
    access the FreezerBox database.  Typically this is done by calling 
    `bind_arguments` on the reactions whenever they're accessed by the app.

    - It is often useful for reactions/reagents to have access to the `appcli` 
      configs used by the app.  This can be done by calling 
      `appcli.share_configs(app, self)` in the `on_bind()` method of argument 
      subclasses.

  - It's important for the structure of the classes to match the structure of 
    the protocol, e.g. a class for each noun.  Because there are so many 
    different interfaces to the app with so much orthogonal functionality, it's 
    really hard to get away with any attempts to fudge this structure.

- The following tests should be included for each protocol:

  - Python API: vary each attribute one at a time, and make sure it has the 
    intended effect.

  - Command-line interface: Provide each command-line argument, and make sure 
    it has the expected effect.  These checks can be much more cursory than for 
    the python API, because we're just making sure that the command-line 
    interface is attached to the protocol correctly; we're already assuming 
    that the protocol itself works.

  - FreezerBox protocol: Check that attributes are correctly inferred from the 
    database, and that multiple FreezerBox reagents can be combined correctly.

    - These test cases are often used twice: once to check that `sw make` 
      works, and once to check that the protocol can be configured from its 
      products.

    - For example, see the PCR and restriction digest tests.

  - FreezerBox attributes: Check that any attributes calculated for the 
    FreezerBox database (e.g. product concentration, product volume, etc.) are 
    correct.

- Note that not all of the scripts in this repository embody all of these 
  recommendations yet.  That's just because it took me a while to settle on the 
  best way to write these scripts, and it takes a lot of work to rewrite all of 
  these scripts.  Any new scripts should adhere to these recommendations, 
  though.
