function [fopt,xopt,gopt]=Polak_Ribiere(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de Polak-Ribière                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de Polak-Ribière";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"0.0005";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------
   compteur = 0;
   compt_vector = [];
   x = xini;
   [compteur, F, G] = Oracle(x, 3, compteur);
   D = - G

   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [compteur, F,G] = Oracle(x,ind, compteur);

//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente
      [compteur, alpha] = Wolfe(alphai, x, D, Oracle, compteur);
      x = x + (alpha*D);
      [compteur, F2, G2] = Oracle(x, ind, compteur);
      
      betak = G2' * (G2 - G) / (norm(G)**2);
      
      D = -G2 + betak * D;

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];
      compt_vector = [compt_vector, compteur]

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de Polak-Ribière')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout, compt_vector);

endfunction
