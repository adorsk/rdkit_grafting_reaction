import unittest

from rdkit import Chem

from rdkit_grafting_reaction import (
    run_grafting_reactions,
    get_graft_site_groups_for_mol,
    make_grafts,
)


class GraftingReactionTestCase(unittest.TestCase):
    def test_run_grafting_reactions(self):
        stock = Chem.MolFromSmarts('[*:1]c[*:2]')
        scion = Chem.MolFromSmarts('[*:1]n[*:2]')
        expected_products_smiles = [
            'c(-n[*:1])[*:1]',
            'c(-n[*:2])[*:2]'
        ]
        products = run_grafting_reactions(stock=stock, scion=scion)
        products_smiles = [Chem.MolToSmiles(product) for product in products]
        self.assertEqual(sorted(products_smiles),
                         sorted(expected_products_smiles))

    def test_get_graft_site_groups_for_mol(self):
        mol = Chem.MolFromSmarts('[*:1]c[*:2][*:1]c[*:3]')
        sites = get_graft_site_groups_for_mol(mol=mol)
        expected_sites = {
            '1': [0, 3],
            '2': [2],
            '3': [5],
        }
        self.assertEqual(sites, expected_sites)

    def test_make_grafts(self):
        stock = Chem.MolFromSmarts('c1cccc1(O[*:1])c[*:1]')
        scion_1 = Chem.MolFromSmarts('c1cccc1([*:1])')
        scion_2 = Chem.MolFromSmarts('O[*:1]')
        product = make_grafts(stock=stock, grafts=[
            {'scion': scion_1, 'stock_idx': 6, 'scion_idx': 5},
            {'scion': scion_2, 'stock_idx': 8, 'scion_idx': 1},
        ])
        product_smiles = Chem.MolToSmiles(product)
        expected_product_smiles = 'Occ1(cccc1)Oc1cccc1'
        self.assertEqual(product_smiles, expected_product_smiles)

    def test_stock_site_groups(self):
        stock = Chem.MolFromSmarts('[*:1][*:2][*:3][*:1]')
        scion = Chem.MolFromSmarts('[*:1][*:4]')
        products = []
        for site_group in [[0], [3], [0, 3]]:
            products.extend(list(run_grafting_reactions(
                stock=stock,
                scion=scion,
                stock_site_groups={'1': site_group}
            )))
        products_smiles = [Chem.MolToSmiles(p) for p in products]
        expected_products_smiles = [
            '[*:1][*:3][*:2][*:4]',
            '[*:1][*:2][*:3][*:4]',
            '[*:2]([*:3][*:4])[*:4]',
        ]
        self.assertEqual(products_smiles, expected_products_smiles)

    def test_scion_site_groups(self):
        stock = Chem.MolFromSmarts('[*:1][*:2][*:3][*:1]')
        scion = Chem.MolFromSmarts('[*:1][*:4][*:5][*:1]')
        products = []
        for site_group in [[0], [3], [0, 3]]:
            products.extend(list(run_grafting_reactions(
                stock=stock,
                scion=scion,
                scion_site_groups={'1': site_group}
            )))
        products_smiles = [Chem.MolToSmiles(p) for p in products]
        expected_products_smiles = [
            '[*:1][*:5][*:4][*:2][*:3][*:4][*:5][*:1]',
            '[*:1][*:4][*:5][*:2][*:3][*:5][*:4][*:1]',
            '[*:1][*:5][*:4][*:2][*:3][*:4][*:5][*:1]',
            '[*:1][*:4][*:5][*:3][*:2][*:4][*:5][*:1]',
            '[*:1][*:4][*:5][*:2][*:3][*:4][*:5][*:1]',
            '[*:1][*:4][*:5][*:2][*:3][*:5][*:4][*:1]',
        ]
        self.assertEqual(products_smiles, expected_products_smiles)

if __name__ == '__main__':
    unittest.main()
