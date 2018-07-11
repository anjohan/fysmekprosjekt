# fysmekprosjekt

# To do:

1. Plot Lennard-Jones potensialkurven. Hva gjør leddene? Kanskje sammenligne med Newtons gravitasjonslov?

2. Finn kraften F ved derivasjon av potensialet, og plot resultatet. Finn ved hvilken avstand r/sigma kraften er 0.

3. Simuler systemet for 2 atomer - skriv xyz-fil, og visualisér med Ovito. Lagre data for alle tidssteg i et 2*3 x nt array, nt = antall tidssteg.
  - Plassér et atom i r1 = [0, 0, 0], og det andre i feks r2 = [1.5, 0, 0] (ved en avstand slik at atomene er fanget i "potensialbrønnen"), og v1 = v2 = [0,0,0]. Hvilken bevegelse forventes, forklar ved hjelp av plottet for potensialet for LJ. Simuler for f.eks. dt = 0.05. Plot både avstand og mekanisk energi som funksjon av tid. I energiplottet vil man sannsynligvis få "spikes". Forklar hvorfor disse forekommer, størrelsen (betydningen) av dem, og hvordan man kan fjerne/redusere disse.
  - Plassér nå atom2 i r2 = [0.95, 0, 0] i stedet. Hva forventes av bevegelse nå? Er energien fortsatt bevart? Begrunn igjen ut fra potensialkurven. Plot de samme plottene.


  -   Implementer forskjellige integrasjonsmetoder, sammenlign plottene for avstand og energi for en gitt dt (ideelt et tidssteg hvor kanskje Euler-Cromer og Vel.Verlet holder oscillasjonene, men hvor Euler sliter). Plot den mekaniske energien som funksjon av tid, kommentér i hvilken grad energien er bevart for de forskjellige metodene, og i hvilke faser av bevegelsen er dette "problematisk" (den mekaniske energien får gjerne noen hopp rundt r = 1 hvor kreftene blir relativt store).
  - Eksperimentér med forskjellige tidssteg for metodene, finn grensene for hver hvor energi fortsatt er bevart.

4. Gjør det samme, men med en metode som kun lagrer to tidssteg om gangen.

5. Mål pot. og kin. energi.

6. Utvid programmet til et vilkårlig antall atomer.
  - Generer fcc gitter
  - Gi atomene tilfeldige initielle hastigheter
  - 2 x 3N arrays for både posisjon og hastighet

  - Gjør flere simuleringer og øk antall atomer gradvis. Rapporter om tidsbruk for simuleringene.
  - Innfør cut-off, og sammenlign igjen tidsbruken

7. Mål energi for f.eks. 512 atomer

8. Implementér forskjellige integrasjonsmetoder (f.eks Euler, Euler-Cromer og Vel.Verlet), og sammenlign disse. Finn maksimalt tidssteg som fortsatt bevarer energi (dette kan også gjøres for 2-atom systemet).

9. Implementer randbetingelser - reflekterende, periodiske

10. Reprodusér resultater fra den sagnomsuste artikkelen fra 60-70 tallet.
