"""
Energy Evaluation Module for RNA Secondary Structures.
Evaluates the free energy of a given dot-bracket structure using ViennaRNA's
thermodynamic model and calculates energy gap metrics.
"""

from typing import Dict
import ViennaRNA as RNA


def evaluate_structure_energy(sequence: str, structure: str) -> float:
    """
    Calculates the Gibbs free energy (kcal/mol) of a given secondary structure
    for a specific RNA sequence using ViennaRNA.

    Args:
        sequence (str): RNA sequence (e.g., "GGAGCAAAAC...").
        structure (str): Dot-bracket structure (e.g., "((((......))))").

    Returns:
        float: Free energy in kcal/mol (rounded to 2 decimal places).
    """
    clean_seq = sequence.upper().replace("T", "U")
    
    # Validation check
    if len(clean_seq) != len(structure):
        raise ValueError(
            f"Sequence length ({len(clean_seq)}) does not match "
            f"structure length ({len(structure)})."
        )

    # Free energy evaluation via ViennaRNA fold_compound
    fc = RNA.fold_compound(clean_seq)
    energy = fc.eval_structure(structure)
    
    return round(energy, 2)


def calculate_energy_gap(predicted_energy: float, reference_mfe_energy: float) -> Dict[str, float]:
    """
    Calculates the absolute energy gap and relative error percentage between 
    a predicted structure (e.g., from VQE) and the classical MFE reference.
    
    Args:
        predicted_energy (float): Energy predicted by the solver (kcal/mol).
        reference_mfe_energy (float): Classical MFE reference energy (kcal/mol).
        
    Returns:
        Dict[str, float]: Energy gap and relative error metrics.
    """
    abs_gap = abs(predicted_energy - reference_mfe_energy)
    
    # Avoid division by zero
    if reference_mfe_energy != 0:
        rel_error = (abs_gap / abs(reference_mfe_energy)) * 100
    else:
        rel_error = 0.0
    
    return {
        "absolute_gap_kcal": round(abs_gap, 2),
        "relative_error_pct": round(rel_error, 2)
    }


if __name__ == "__main__":
    # Test on Moderna official 44 nt sequence
    official_seq = "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG"
    mfe_struct   = "((((((.((((((((...))))))))(((((...))))))))))"
    mfe_energy   = -13.10
    
    # Candidate structure test (e.g., slightly sub-optimal quantum prediction)
    test_struct  = "((((((.((((((......))))))(((((...))))))))))"
    test_energy  = evaluate_structure_energy(official_seq, test_struct)
    
    gap_metrics  = calculate_energy_gap(test_energy, mfe_energy)

    print("=== ENERGY EVALUATION BENCHMARK TEST ===")
    print(f"Sequence         : {official_seq[:15]}...")
    print(f"Candidate Struct : {test_struct}")
    print(f"Predicted Energy : {test_energy} kcal/mol")
    print(f"Reference MFE    : {mfe_energy} kcal/mol")
    print(f"Absolute Gap     : {gap_metrics['absolute_gap_kcal']} kcal/mol")
    print(f"Relative Error   : {gap_metrics['relative_error_pct']} %")
