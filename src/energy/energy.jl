function energy_naive(interaction::FSSoG_naive{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position.coo for i in 1:interaction.n_atoms]

    E_long = long_energy_naive(interaction, position, charge)
    E_short = short_energy_naive(interaction, neighbor, position, charge)

    return E_long + E_short
end

function ExTinyMD.energy(interaction::FSSoGInteraction{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    for i in 1:interaction.n_atoms
        interaction.charge[i] = atoms[info.particle_info[i].id].charge
        interaction.position[i] = info.particle_info[i].position.coo
    end

    E_long = energy_long(interaction)
    E_mid = energy_mid(interaction)
    E_short = energy_short(interaction, neighbor)

    return E_long + E_mid + E_short
end