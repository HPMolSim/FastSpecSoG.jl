function energy_naive(interaction::FSSoG_naive{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:interaction.n_atoms]

    energy_long = long_energy_naive(interaction, position, charge)
    energy_short = short_energy_naive(interaction, neighbor, position, charge)

    return energy_long + energy_short
end