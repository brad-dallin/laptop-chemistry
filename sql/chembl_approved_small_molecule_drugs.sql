SELECT DISTINCT m.chembl_id AS chembl_molecule_id,
s.canonical_smiles AS smiles,
-- Molecule Dictionary
m.molregno, m.pref_name, m.max_phase, m.therapeutic_flag, m.dosed_ingredient,
m.structure_type, m.chebi_par_id, m.molecule_type, m.first_approval, m.oral,
m.parenteral, m.topical, m.black_box_warning, m.natural_product, m.first_in_class,
m.chirality, m.prodrug, m.inorganic_flag, m.usan_year, m.availability_type,
m.usan_stem, m.polymer_flag, m.usan_substem, m.usan_stem_definition,
m.indication_class, m.withdrawn_flag, m.chemical_probe, m.orphan,
-- Compound Properties
p.mw_freebase, p.alogp, p.hba, p.hbd, p.psa, p.rtb, p.ro3_pass,
p.num_ro5_violations, p.cx_most_apka, p.cx_most_bpka, p.cx_logp,
p.cx_logd, p.molecular_species, p.full_mwt, p.aromatic_rings,
p.heavy_atoms, p.qed_weighted, p.mw_monoisotopic, p.full_molformula,
p.hba_lipinski, p.hbd_lipinski, p.num_lipinski_ro5_violations,
p.np_likeness_score
-- Join tables
FROM compound_structures s
    RIGHT JOIN compound_properties p ON s.molregno = p.molregno
    RIGHT JOIN molecule_dictionary m ON p.molregno = m.molregno
    JOIN compound_records r ON m.molregno = r.molregno
    AND m.max_phase = 4
    AND m.molecule_type = 'Small molecule'
    AND m.inorganic_flag = 0
ORDER BY m.first_approval ASC