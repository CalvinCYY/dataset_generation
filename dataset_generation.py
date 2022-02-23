from rdkit import Chem
from rdkit.Chem import AllChem
import pybel as pyb

from tqdm import tqdm

def batch_process_sdf(filepath, outfolder, dataset_name="", sdf_id_tag=None, make_3D=True, progress=True):
    suppl = Chem.SDMolSupplier(filepath)
    if progress:
        suppl = tqdm(suppl)

    mol_list = []
    for i, rdmol in enumerate(suppl):
        if not rdmol.HasProp('_Name'):
            if not sdf_id_tag:
                rdmol.SetProp('_Name', f'{dataset_name}{i}')
            else:
                rdmol.SetProp('_Name', sdf_id_tag)

        if make_3D:
            processed_mol = make_3D_rdkit(rdmol)

        else:
            processed_mol = rdmol

        molid = processed_mol.GetProp('_Name')

        file_id = f'{dataset_name}_{molid}'

        writer = Chem.SDWriter(f'{outfolder}/{file_id}.sdf')
        writer.write(processed_mol)
        writer.flush()

def make_3D_rdkit(rdmol):
    rdmol_3d = Chem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol_3d)
    return rdmol_3d