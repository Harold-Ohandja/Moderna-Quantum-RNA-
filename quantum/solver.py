from typing import Dict, List, Tuple

def brute_force_qubo_minimum(
    qubo: Dict[Tuple[int, int], float],
    num_vars: int
) -> Tuple[List[int], float]:
    """Find the minimum QUBO cost by brute force for tiny instances."""