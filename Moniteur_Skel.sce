///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_R.sce');
   exec('Structures_R.sce');

   // Affichage des resultats

   exec('Visualg.sci');

   // Verification  des resultats

   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   exec('Oracle.sci');
   exec('Oracle_lagrange.sci')
   exec('Optim_Scilab.sci');
   exec('Wolfe.sci');
   exec('Gradient_W.sci');
   exec('Gradient_F.sci');
   exec('Newton.sce')
   exec('Polak_Ribiere.sci')
   exec('BFGS.sci')
   //titrgr = "Fonction optim de Scilab sur le probleme primal";

   // -----> A completer...
   // -----> A completer...
   // -----> A completer...

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal

//   xini = 0.1 * rand(n-md,1);
   xini = 0.1 * rand(md,1);

// ----------------------------0
// Minimisation proprement dite
// ----------------------------

   // Exemple : la fonction "optim" de Scilab
   //
//   [fopt,xopt,gopt] = Gradient_F(OraclePH,xini);
//   [fopt,xopt,gopt] = Gradient_W(OraclePH,xini);
//    [fopt,xopt,gopt] = Polak_Ribiere(OraclePH,xini);
//     [fopt,xopt,gopt] = BFGS(OraclePH,xini);
   // -----> A completer...
//         [fopt,xopt,gopt] = Gradient_F(OraclePHL,xini);
//   [fopt,xopt,gopt] = Gradient_W(OraclePHL,xini);
//   [fopt,xopt,gopt] = Polak_Ribiere(OraclePHL,xini); // 1. 408
//     [fopt,xopt,gopt] = BFGS(OraclePHL,xini); // 0.22
//    [fopt,xopt,gopt] = Newton(OraclePHL,xini); // 0.052
// --------------------------
// Verification des resultats
// --------------------------

//   [q,z,f,p] = HydrauliqueP(xopt);

   [q,z,f,p] = HydrauliqueD(xopt);

   Verification(q,z,f,p);

//
