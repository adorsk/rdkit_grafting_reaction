# rdkit_grafting_reaction
A library for joining molecular fragments together with RDKit.

ACHTUNG! This library is still new, and has not been widely tested. Caveat programmer.

## Analogy: Grafting
This library uses the analogy of [grafting](https://en.wikipedia.org/wiki/Grafting).

In a horticultural graft a branch from one plant is joined to another plant. The branch being donated is called the 'scion', and the receiving plant is called the 'rootstock' or the 'stock'.

A stock can be grafted with multiple scions, at difference locations.

In this library we use this analogy to define reactions for combine molecular fragments.

One molecule acts as the stock. Multiple scion fragments can be joined to the stock molecule.

## run_grafting_reactions
The primary method in this library is the `run_grafting_reactions` method.

A grafting reaction takes four parameters:

1. A stock molecule.
2. A scion molecule.
3. (optional) A set of grafting site groups on the stock molecule.
4. (optional) A set of grafting sites groups on the scion molecule.

If grafting site groups are not specified, we use all available grafting sites to make combinations of grafts.

### Grafting Sites
Molecules can have multiple groups of grafting sites. A grafting site is a point where one molecule can be joined to another. This library uses SMILES atom maps to define grafting sites.

For example, for the molecule with SMILES `c[*:1]`, the grafting site is the `[*:1]` part.

Molecules can have groups of grafting sites. For example, the SMILES `c[*:1][*:2][*:1]` contains two groups of grafting sites: the group with atom map label `1`, and the group with atom map label `2`. The group for `1` contains two grafting sites. The purpose of grafting site groups is to constrain where certain fragments can attach.

For example, you have may have a library of bridge fragments, with two possible grafting sites. On one site, you want to graft a core fragment, and on the other, you want to graft a terminal fragment. You could use different grafting site groups to limit graft sites to specific fragment type.
