function []=Visualg(logG,logP,Cout, compt_vector)


///////////////////////////////////////////////////////////////
//                                                           //
//   VISUALISATION DES ITEREES D'UN ALGORITHME               //
//                                                           //
// Affichage des resultats obtenus durant l'algorithme       //
//                                                           //
//   - logG : tableau contenant le log base 10 de la norme   //
//            du gradient du critere pour chaque iteration   //
//            de l'algorithme d'optimisation.                //
//   - logP : tableau contenant le log base 10 du pas de     //
//            gradient pour chaque iteration.                //
//   - Cout : tableau contenant la valeur du critere pour    //
//            chaque iteration.                              //
//                                                           //
// L'affichage du cout est probablement le moins interessant //
// et peut etre omis. Afficher le log en base 10 de la norme //
// du gradient permet bien representer la convergence vers 0 //
// du gradient. Enfin,  afficher le log en base 10 du pas de //
// descente donne par l'algorithme de recherche lineaire est //
// un bon  moyen de verifier que cet algorithme de recherche //
// lineaire fonctionne de maniere satisfaisante.             //
//                                                           //
///////////////////////////////////////////////////////////////


numwin = 10;
typvis =  1;

[nlig,ncol] = size(logG);
absXpas = [1:nlig]';
absX = [1:nlig]';
absX = compt_vector';
disp(size(compt_vector'))
disp(size(logG))

if typvis == 0 then

//
// Affichage des 3 graphiques, un par fenetre
//

   xset("window",numwin);
   clf(numwin);
   plot2d(logG);
   xtitle('Norme du gradient (echelle logarithmique)','Iter.','log||G||');

   numwin = numwin + 1;
   xset("window",numwin);
   clf(numwin);
   plot2d(absXpas, logP);
   xtitle('Pas de gradient (echelle logarithmique)','Iter.','log(alpha)');

   numwin = numwin + 1;
   xset("window",numwin);
   clf(numwin);
   plot2d(Cout);
   xtitle('Evolution du critere','Iter.','Cout');

else

//
// Affichage de 2 graphiques sur une meme fenetre
//

   xset("window",numwin);
   xset("wdim",600,970);
   clf(numwin);

//   xtitle(titrgr);

  subplot(311);
  plot2d(absXpas,logG,style=5,...
        leg='Norme du gradient (echelle logarithmique) '+...
           'en fonction des itérations');

   subplot(312);
   plot2d(absX,logG,style=5,...
          leg='Norme du gradient (echelle logarithmique) '+...
              'en fonction des appels à l oracle');

   subplot(313);
   plot2d(absXpas,logP,style=2);

end

endfunction
