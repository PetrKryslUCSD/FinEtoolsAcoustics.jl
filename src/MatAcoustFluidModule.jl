"""
    MatAcoustFluidModule

Module for acoustic-fluid  material.
"""
module MatAcoustFluidModule

import FinEtools.MatModule: AbstractMat

"""
    MatAcoustFluid <: AbstractMat

Type for acoustic fluid material.
"""
struct MatAcoustFluid{T} <: AbstractMat
	bulk_modulus::T;# Bulk modulus
	mass_density::T;# Mass density
end

"""
    bulkmodulus(self::MatAcoustFluid)

Return the bulk modulus.
"""
function bulkmodulus(self::MatAcoustFluid)
	return self.bulk_modulus
end

end
