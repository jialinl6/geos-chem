import DocStringExtensions

"""
Contains derived types and methods to create a registry of each variable contained within a given module.  This will allow the user to obtain a pointer to each module variable by searching for its name.
"""
module registry_mod

export RegItem

"""
Derived type for a REGISTRY ITEM (a single registry entry). This represents a single module variable, plus some metadata.
"""
struct RegItem
	# Identifying info
	"e.g. \"STATE_VARIABLE\""
  full_name::String
	"Name of state"
  state::String
	"Name of variable"
  variable::String

	# Metadata
	"Longer description"
	description::String
	"Memory use in Kb"
	memory_in_kb::AbstractFloat
	"Numerical KIND value of data"
	source_kind_val::Integer
	"Numerical KIND value for output"
	output_kind_val::Integer
	"Dimensions of data"
	rank::Integer
	"Units of data"
	units::String
	"e.g. \"xyz\", \"yz\", \"y\", \"t\""
	dim_names::String
	"Is data on level edges (T/F)?"
	on_level_edges::Bool

	# Pointers to floating point data (8-byte precision)
	"For 0D 8-byte data"
  ptr_0d_8::AbstractFloat
	"For 1D 8-byte data"
  ptr_1d_8::Vector{AbstractFloat}
	"For 2D 8-byte data"
  ptr_2d_8::Matrix{AbstractFloat}
	"For 3D 8-byte data"
  ptr_3d_8::Array{AbstractFloat, 3}

	# Pointers to floating point data (4-byte precision)
	"For 0D 4-byte data"
	ptr_0d_4::AbstractFloat
	"For 1D 4-byte data"
	ptr_1d_4::Vector{AbstractFloat}
	"For 2D 4-byte data"
	ptr_2d_4::Matrix{AbstractFloat}
	"For 3D 4-byte data"
	ptr_3d_4::Array{AbstractFloat, 3}

	# Pointers to integer data
	"For 0D int data"
	ptr_0d_I::Integer
	"For 1D int data"
	ptr_1d_I::Vector{Integer}
	"For 2D int data"
	ptr_2d_I::Matrix{Integer}
	"For 3D int data"
	ptr_3d_I::Array{Integer, 3}
end

"""
Derived type for a METAREGISTRY ITEM (a linked-list of REGISTRY ITEMS)
"""
struct MetaRegItem
	"Pointer to next node"
	next::MetaRegItem
	"Registry item within"
	item::RegItem
end

end
