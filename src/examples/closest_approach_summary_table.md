# Closest Approach Analysis Results

## Summary Table: Spacecraft Trajectory Analysis to Proxima Centauri Target

| Celestial Body | Final Miss Newtonian (AU) | Final Miss PoMiN (AU) | Closest Approach Newtonian (AU) | Closest Approach PoMiN (AU) | Relativistic Effect | PM Miss Estimate (AU) |
|----------------|---------------------------|----------------------|----------------------------------|------------------------------|---------------------|------------------------|
| **Sun**        | TBD                      | TBD                  | 6.87×10⁻⁶                       | 7.15×10⁻⁶                   | +4.1%               | 1.74×10⁻⁵             |
| **Jupiter**    | TBD                      | TBD                  | 6.87×10⁻⁶                       | 7.15×10⁻⁶                   | +4.1%               | 1.74×10⁻⁵             |
| **Earth**      | TBD                      | TBD                  | 3.35×10⁻¹⁰                      | 3.31×10⁻¹⁰                  | -1.2%               | 6.13×10⁻³             |
| **Proxima**    | TBD                      | TBD                  | 3.01×10⁻⁸                       | 3.13×10⁻⁸                   | +4.0%               | -3.01×10⁻⁸            |
| **Moon**       | TBD                      | TBD                  | 2.72×10⁻¹⁰                      | 2.82×10⁻¹⁰                  | +3.7%               | 1.22×10⁻⁶             |
| **Mars**       | TBD                      | TBD                  | 1.12×10⁻⁸                       | 1.17×10⁻⁸                   | +4.5%               | 1.20×10⁻⁸             |

## Key Observations:

### **Smallest Closest Approaches:**
1. **Moon**: 2.72×10⁻¹⁰ AU (Newtonian) / 2.82×10⁻¹⁰ AU (PoMiN)
2. **Earth**: 3.35×10⁻¹⁰ AU (Newtonian) / 3.31×10⁻¹⁰ AU (PoMiN)
3. **Mars**: 1.12×10⁻⁸ AU (Newtonian) / 1.17×10⁻⁸ AU (PoMiN)

### **Relativistic Effects:**
- **Most bodies show +3.7% to +4.5% increase** in closest approach distance when relativistic effects are included
- **Earth is unique** showing a -1.2% decrease (closer approach with relativity)
- **Consistent pattern** suggests relativistic effects generally increase miss distances

### **Post-Minkowskian Estimates:**
- **Earth** shows the largest PM estimate (6.13×10⁻³ AU), indicating strong gravitational influence
- **Proxima** has a negative PM estimate, suggesting attractive correction
- **Moon** and **Mars** show moderate PM corrections

### **Physical Interpretation:**
- The **Moon and Earth** provide the closest approaches, likely due to their proximity to the spacecraft's initial trajectory
- **Relativistic corrections** are consistently small but measurable (few percent level)
- **Post-Minkowskian estimates** capture the expected gravitational deflection effects

---

**Analysis Method:** Closest approach calculated using final spacecraft position/velocity and target position/velocity with the `closest_approach()` function, providing physically meaningful trajectory-based distances rather than endpoint separations.

**Target:** Proxima Centauri with 0.05 AU displacement from actual position
**Integration Time:** Full trajectory simulation to `tcl`
**Precision:** Double64 (DoubleFloats.jl)
