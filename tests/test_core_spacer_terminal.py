import itertools
import os
import unittest

from rdkit import Chem

from rdkit_grafting_reaction import (
    get_graft_site_groups_for_mol,
    run_grafting_reactions,
)

THIS_DIR = os.path.dirname(__file__)

def squash_duplicate_mols(mols):
    return [Chem.MolFromSmiles(smiles)
            for smiles in set([Chem.MolToSmiles(mol) for mol in mols])]

class CoreSpacerTerminalTestCase(unittest.TestCase):
    def test_acceptors(self):
        fragments_smiles = {
            'cores': [
                '[He]C1=C([He])C2=CC(C#N)=CC=C2S1',
                'O=C1C2=C([He])SC([He])=C2C(C3=C(C(F)(F)F)SC(C(F)(F)F)=C31)=O',
            ],
            'spacers': [
                '[He]C1=CN=C([He])N=N1',
                '[He]C1=C2C(OCCO2)=C([He])[Si](C)1C',
            ],
            'terminals': [
                '[He]C#CC1=C(F)C(F)=C(C#N)S1',
                'O=C(O[He])C1=C(C#N)ON=C1C#N',
            ],
        }

        def replace_handles_with_atom_maps(smiles):
            handles = ['[He]']
            for i, handle in enumerate(handles):
                smiles = smiles.replace(handle, '[*:%s]' % (i + 1))
            return smiles

        fragments_mols = {
            key: [Chem.MolFromSmiles(replace_handles_with_atom_maps(smiles))
                  for smiles in smiles_list]
            for key, smiles_list in fragments_smiles.items()
        }

        # Graft spacers on to cores.
        cores_spacers = []
        for core, spacer in itertools.product(fragments_mols['cores'],
                                              fragments_mols['spacers']):
            # Leave one spacer site free for the terminal.
            spacer_site_groups = get_graft_site_groups_for_mol(spacer)
            spacer_sites = spacer_site_groups['1']
            for spacer_site_combos in itertools.combinations(
                spacer_sites, len(spacer_sites) - 1):
                for core_spacer in run_grafting_reactions(
                    stock=core,
                    scion=spacer,
                    scion_site_groups={'1': spacer_site_combos}
                ):
                    cores_spacers.append(core_spacer)
        cores_spacers = squash_duplicate_mols(cores_spacers)

        # Graft terminals onto (core+spacer) combos.
        cores_spacers_terminals = []
        for terminal, core_spacer in itertools.product(
            fragments_mols['terminals'], cores_spacers):
            for core_spacer_terminal in run_grafting_reactions(
                stock=core_spacer,
                scion=terminal
            ):
                cores_spacers_terminals.append(core_spacer_terminal)
        cores_spacers_terminals = squash_duplicate_mols(cores_spacers_terminals)

        actual_smiles = [Chem.MolToSmiles(core_spacer_terminal)
                         for core_spacer_terminal in cores_spacers_terminals]
        expected_smiles_path = os.path.join(
            THIS_DIR, 'expected_core_spacer_terminal_smiles.txt')
        with open(expected_smiles_path) as f:
            expected_smiles = [l.strip() for l in f.readlines()]
        self.assertEqual(sorted(actual_smiles), sorted(expected_smiles))
