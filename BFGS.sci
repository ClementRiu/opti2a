function [fopt,xopt,gopt]=BFGS(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de BFGS                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de BFGS";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"0.0005";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   compteur = 0;
   compt_vector = [];
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   W = eye(x * x')

   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [compteur, F1, G1] = Oracle(x,ind, compteur);

//    - test de convergence

      if norm(G1) <= tol then
         kstar = k;
         break
      end
      

//    - calcul de la direction de descente

      D = - W * G1;

//    - calcul de la longueur du pas de gradient

      [compteur, alpha] = Wolfe(alphai, x, D, Oracle, compteur);

//    - mise a jour des variables
      x1 = x;
      x = x + (alpha * D);
      [compteur, F2, G2] = Oracle(x, ind, compteur)

      deltax = x -x1;
      deltag = G2 - G1
      deltaxg = deltax * deltag';
      deltagx = deltag * deltax';
      deltaxx = deltax * deltax';
      deltagx_scalar = deltax' * deltag;
      W = (eye(deltaxg) - deltaxg / deltagx_scalar) * W * (eye(deltagx) - deltagx / deltagx_scalar) + deltaxx / deltagx_scalar;
//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G1)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F1 ];
      compt_vector = [compt_vector, compteur];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F1;
   xopt = x;
   gopt = G2;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode BFGS')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout, compt_vector);

endfunction
