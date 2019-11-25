.. AntimonyCombinations documentation master file, created by
   sphinx-quickstart on Sun Nov 24 10:57:34 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AntimonyCombinations's documentation!
================================================

`AntimonyCombinations` is a package developed on top of
[tellurium](http://tellurium.analogmachine.org/) and
[antimony](http://antimony.sourceforge.net/) for building
[sbml](http://sbml.org/Main_Page) models in a combinatorial
way.

The idea is that you have a core model which you
are more confident in regarding its structure and an arbitrary
number of additional hypotheses, called hypothesis extensions.
`AntimonyCombinations` provides a way of quickly building the
comprehensive set of model topologies, given the core hypothesis
and hypothesis extensions.

Installation
------------

    `$ pip install AntimonyCombinations`


Import
------

The pip command uses camel caps but you still import the code using the
Python convention of separation by `'_'` for package names.

.. code-block::
   :linenos:

   from antimony_combinations import Combinations, HypothesisExtension

.. toctree::
   :maxdepth: 2

   combinations.rst
   hypothesis_extension.rst

