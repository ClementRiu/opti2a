function [compteur, alphan,ok]=Wolfe(alpha,x,D,Oracle, compteur)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1;
   omega2 = 0.9;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [compteur, F,G] = Oracle(x,ind, compteur)

   // Initialisation de l'algorithme

   alphan = alpha;
   xn     = x;
   xp     = 2 * abs(x) + abs(dltx) 

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

   while ok == 0

      xp = xn;
      xn = x + (alphan*D);

      // Calcul des conditions de Wolfe

      [compteur, F2,G2] = Oracle(x,4, compteur);
      [compteur, F3,G3] = Oracle(xn,4, compteur);


      cond1 = F3 <= F2 + omega1 * alphan * G2' * D;
      cond2 = G3' * D >= omega2 * G2' * D;

      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      //   faire ok = 1 : on sort alors de la boucle while
      // - sinon, modifier la valeur de alphan : on reboucle.

      
      if ~cond1 then
         alphamax = alphan;
         alphan =  (alphamin + alphamax) / 2;
      else
         if ~cond2 then
            alphamin = alphan;
            if alphamax == %inf then
               alphan = 2 * alphamin;
            else
               alphan =  (alphamin + alphamax) / 2;
            end
         else 
            ok = 1;
         end
      end

      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end

endfunction
