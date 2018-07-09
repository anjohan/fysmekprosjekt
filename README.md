# fysmekprosjekt

# To do:

1. Plot Lennard-Jones potensialkurven. Hva gjør leddene?

2. Finn kraften F ved derivasjon av potensialet, og plot resultatet.

3. Simuler systemet for 2 atomer - skriv xyz-fil, og visualisér med Ovito. Lagre data for alle tidssteg i et 2*3 x nt array, nt = antall tidssteg.

4. Gjør det samme, men med en metode som kun lagrer to tidssteg om gangen.

5. Mål pot. og kin. energi.

6. Utvid programmet til et vilkårlig antall atomer.
  - Generer fcc gitter
  - Gi atomene tilfeldige initielle hastigheter
  - 2 x 3N arrays for både posisjon og hastighet

7. Mål energi for f.eks. 512 atomer

8. Implementér forskjellige integrasjonsmetoder (f.eks Euler, Euler-Cromer og Vel.Verlet), og sammenlign disse. Finn maksimalt tidssteg som fortsatt bevarer energi (dette kan også gjøres for 2-atom systemet).

9. Implementer randbetingelser - reflekterende, periodiske

10. Reprodusér resultater fra den sagnomsuste artikkelen fra 60-70 tallet.
