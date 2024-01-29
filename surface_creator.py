# Import packages
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
import numpy as np

def surface_generator(thickness, width, depth, rhombohedral_miller, padding=1.0):
    """This function generates an orhorhombic hematite Fe2O3 layer with 
    specified Miller plane exposed and of desired dimensions and padding.
    
    Parameters
    ------------
        thickness: float
            The minimum thickness of the hematite Fe2O3 layer in Angstroms.
        width: float
            The required width of the hematite Fe2O3 layer in Angstroms.
        depth: float
            The required depth of the hematite Fe2O3 layer in Angstroms.
        rhombohedral_miller: tuple
            The Miller indices of the exposed plane. Should be format (h,k,l)
        padding: float, Optional
            The vacuum on either side of the layer.
    """

    # Create Structure object of orthorhombic Fe2O3 structure.
    fe2o3 = Structure.from_file("Fe2O3_orth.cif")

    # BRIEF ASIDE ON MILLER THE NEED FOR MILLER INDEX CONVERSION:
    #
    #   1. The "Fe2O3_orth.cif" file contains an orthorhombic cell 
    #       made by taking the rhomohedral unit cell of hematite 
    #       Fe2O3, creating a 1x2x1 supercell and transforming 
    #       it onto an orthorhombic lattice.
    #   2. This orthorhombic unit cell makes creating 
    #       orthorhombic layers easier.
    #   3. Miller indices, however, must be converted from the 
    #       rhombohedral system to the orthorhombic one.
    #
    #   100:
    #       In a rhombohedral system, (100), (010), (1-10) 
    #       are equivalent. When the system is converted to 
    #       orthorombic, by extending in the "b" direction, 
    #       these correspond to (110), (010) and (1-10). 
    #       In this conversion, all are converted to (110).
    #
    #   001:
    #       The (001) plane is the same in 
    #       both rhombohedral and orhtorombic systems.
    #
    #   110:
    #       In a rhombohedral system, (110), (1-20), (2-10)
    #       are equivalent. When the system is converted to
    #       orthorhombic, by extending in the "b" direction,
    #       these correspond to the (100), (130), (1-30).
    #       In this conversion, all are converted to (100).
    #
    #   012:
    #       In a rhombohedral system, (012), (-102), (1-12)
    #       are equivalent. When the system is converted to
    #       orthorhombic, by extending in the "b" direction,
    #       these correspond to the (022), (1-12), (-1-12).
    #       In this conversion, all are converted to (022).

    # Convert rhombohedral to orthorombic Miller indices.
    plane_100 = [(1,0,0), (0,1,0), (1,-1,0), (-1,0,0), (0,-1,0), (-1,1,0)]
    plane_001 = [(0,0,1), (0,0,-1)]
    plane_110 = [(1,1,0), (1,-2,0), (2,-1,0), (-1,-1,0), (-1,2,0), (-2,1,0)]
    plane_012 = [(0,1,2), (1,-1,2), (-1,-1,2), (0,-2,-2), (-1,1,-2), (1,1,-2)]

    if rhombohedral_miller in plane_100:
        mil = (1,1,0)
    elif rhombohedral_miller in plane_001:
        mil = (0,0,1)
    elif rhombohedral_miller in plane_110:
        mil = (1,0,0)
    elif rhombohedral_miller in plane_012:
        mil = (0,2,2)
    
    # Generate layers.
    slabgen = SlabGenerator(fe2o3,
                    miller_index=mil,
                    min_slab_size=thickness,
                    min_vacuum_size=padding,
                    center_slab=True,
                    lll_reduce=True,
                    primitive=False)
    
    slabs = slabgen.get_slabs()

    # Scale up to required dimensions and export.
    mult_a = int(round(width/slabs[0].lattice.a, 0))
    mult_b = int(round(depth/slabs[0].lattice.b, 0))

    for s, i in zip(slabs, np.arange(0, len(slabs))):
        s.make_supercell([mult_a, mult_b, 1])
        s.get_orthogonal_c_slab()
        s.to(filename=f'slab_{rhombohedral_miller}_{i}.cif')
