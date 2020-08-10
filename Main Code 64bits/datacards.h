// PERMET DE LIRE DES DATA A PARTIR D'UN FICHIER.
// Ecrit par:
// Eric Aubourg, Decembre 95 // DAPNIA/SPP (Saclay) / CEA    LAL - IN2P3/CNRS  (Orsay)
// Reza Ansari, Aout 96
// Daniel Comparat Sept 2012


//   	Cette classe permet la gestion des parametres d'un programme a partir
//   	de mot-cle (lecture d'un fichier par exemple)


// Les "card" sont les variables, ce sont les noms en début de ligne commençant par @  suivit ensuite une liste de valeurs séparées par des blancs (ou tab) et ceux jusqu'a un saut de ligne

/*** Example de fichier datacard **********

#
@SEEING_PRCTAGE 0.9
@TRUC_ENTIER 15 17 23
@TRUC_STRING udchdjkzbcz
#

*************************/

// La 1ère valeur peut être int (appelée par: data.IParam("nom_de_la_variable")) ou double (de même avec D) où string (S);
// De même si on veux la troisième valeur data.IParam("nom_de_la_variable",2)
// Toute autre mot ou ligne est ignorée


// int   NbCards()
//	Renvoie le nombre de cartes data
// bool	 HasKey(string const& key)
//	Indique l'existence d'une carte avec la cle "key"
// int   NbParam(string const& key)
//	Indique le nombre de parametres (separes par des espaces) pour la cle "key"



/******************  Exemple Utilisation datacard ******************


string nomdat = "ma_data_card" ;
DataCards data(nomdat.c_str());
double prctage =  data.DParam("SEEING_PRCTAGE") ;
int itruc = data.IParam("TRUC_ENTIER") ;
string struc = data.SParam("TRUC_STRING") ;


***************************/

#ifndef DATACARDS_SEEN
#define DATACARDS_SEEN

#include <string>
#include <functional>
#include <list>
#include <vector>
#include <string>

using namespace std;


typedef int (*ProcCard)(string const& key, string const& toks);

class DataCards  {
public:
  DataCards();
  DataCards(string const& fn);

  virtual ~DataCards() {}

// AddProcF(ProcCard f, string const& mtch="*")
//	Ajoute une fonction de traitement a la liste pour les mots cle
//      compatibles avec la chaine "mtch" ("mtch" peut contenir "*" en debut
//      fin de mot)
//


//	Rajoute la carte "line" a la liste
  void    AddProcF(ProcCard f, string const& mtch="*");


//	Supprime les cartes existantes
  void    Clear();

//      Lit le contenu du fichiers "fn" et ajoute les cartes a la liste
  void    ReadFile(string const& fn);
  void    AppendCard(string const& line);

  int     NbCards();
  bool    HasKey(string const& key);
  int     NbParam(string const& key);
  string  SParam(string const& key, int numero = 0, string def="");
  long    IParam(string const& key, int numero = 0, long def = 0);
  double  DParam(string const& key, int numero = 0, double def = 0);

  friend ostream& operator << (ostream& s, DataCards c);

public:
  struct Card {
    string kw;
    vector<string> tokens;
    //	 STRUCTCOMPF(Card,kw)
  };
  typedef list<Card> CardList;
  struct CrdPF {
    ProcCard pf;
    string  patt;
    // STRUCTCOMPF(CrdPF,pf)
  };
  typedef list<CrdPF> CrdPFList;
public:
  CardList cards;
  CrdPFList cpfs;

  void  DoReadFile(string const& fn);

  int   ApplyPF(CrdPF & cpf, string const& key, string const& toks);
  int   ApplyPFL(string const& key, string const& toks);

  void  RemoveCard(string const& key);

  Card* FindKey(string const& key);

#ifndef SWIG
 struct KeyEq : binary_function<Card, string, bool> {
   bool operator()(const Card& x, const string& y) const
   { return x.kw == y; }
  };
#endif

};
#endif
