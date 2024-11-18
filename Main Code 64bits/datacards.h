
// ENABLES READING DATA FROM A FILE.
// Written by:
// Eric Aubourg, December 1995 // DAPNIA/SPP (Saclay) / CEA LAL - IN2P3/CNRS (Orsay)
// Reza Ansari, August 1996
// Daniel Comparat, September 2012

// This class manages program parameters using keywords (e.g., reading from a file).

// "Cards" refer to variables, which are the names at the beginning of a line starting with @,
// followed by a list of values separated by spaces (or tabs) until the end of the line.


/*** Example de fichier datacard **********

#
@SEEING_PRCTAGE 0.9
@TRUC_ENTIER 15 17 23
@TRUC_STRING udchdjkzbcz
#

*************************/

// The first value can be an integer (accessed via: `data.IParam("variable_name")`),
// a double (accessed via `data.DParam`), or a string (`data.SParam`).
// To access the third value, for example: `data.IParam("variable_name", 2)`.
// Any other word or line is ignored.

// int NbCards()
//    Returns the number of data cards.
// bool HasKey(string const& key)
//    Indicates whether a card with the given key exists.
// int NbParam(string const& key)
//


/******************  Exemple Utilisation datacard ******************


string datacard_file = "my_data_card";
DataCards data(datacard_file.c_str());
double percentage = data.DParam("SEEING_PRCTAGE");
int integer_value = data.IParam("TRUC_ENTIER");
string string_value = data.SParam("TRUC_STRING");


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
