/** Classe Param�tre (Param) et Liste de Param�tre (FitParams)

Un param�tre est:
un nom
un bool�an is_scanned disant si il est fix� ou non. Par d�faut false
un bool�an is_time_dependent disant si il est varie dans le temps ou pas. Par d�faut false
Une liste de:
 val,  minv, maxv, nbstep, tau, val_t0;
 // val_t0 = minv + steps * (maxv - minv )/nbstep
 // tau =   taux de variation temporelle --> val = val_t0 exp^(-t/tau)


Exemple d'utilisation de FitParam


  FitParams &params; // this is a vector<Param *>
  Param *pDayMax;
  pDayMax  = params.LocateOrAddNonConstParam("DayMax");


// fixer un parametre  de nom paramname
  param_a_fixer->is_scanned = false;

// mettre un parametre de nom paramname a une certaine valeur value
  Param* param_to_set =  params.LocateOrAddParam(paramname);
    param_to_set->val = value; //peut rajouter  pour v�rifier si il est fix�: if(!param_to_set->is_scanned)


// pour initialiser tout d'un coup avec un tableau tab de la bonne taille :
 for(size_t i=0; i<params.size(); ++i) params[i]->val=tab[i];



// boucle sur les parametres de params :
  for (ParamIterator i=params.begin(); i != params.end(); ++i)
    {
      Param &p = **i;
	if (p.is_scanned) continue;
	// faire la boucle sur la valeur de param ici
	for (int n = 0 ; n < p.nstep ; n++)
	// etc.

      }

****************************/

#ifndef PARAM__H
#define PARAM__H

#include <iostream>
#include <vector>
#include <iostream>
#include <string>

#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for flat random generator

using namespace std;

//! Named parameters (i.e. what estimation theory call parameters).
class Param
{
public :

    double val,  minv, maxv,  val_t0, tau;
    int nbstep;
    bool is_scanned, is_time_dependent;
    std::string name;
    Param() : val(0.), minv(0.), maxv (0.), val_t0(0.), tau(1e100), nbstep(1), is_scanned(false), is_time_dependent(false), name("UNKNOWN") {};
    Param(const string &Name) : val(0.),  minv(0.), maxv (0.), val_t0(0.),  tau(1e100), nbstep(1), is_scanned(false), is_time_dependent(false), name(Name) {};

//  Donne le nom.
    std::string Name() const
    {
        return name;
    }

// Sortie du param�tre (val, minv , ...)
    void dump(ostream &stream = cout, bool nice_formated_output=false) const;

    friend ostream& operator << (ostream &stream, const Param &p)
    {
        p.dump(stream);
        return stream;
    }
};


//! the ensemble of all parameters of the fit
// It is vector of pointers, in order to keep always the proper allocation.

class FitParams : public vector<Param *>
{
    string name;

public :

    FitParams(const string &Name) : name(Name) {};
    FitParams() : name("NoName") {};
    void SetName(const string &Name)
    {
        name=Name;
    }

    FitParams(const FitParams &);

//pour atteindre un param Retourne le pointeur Param* vers le param qui � le nom ParamName
    const Param* LocateParam(const string &ParamName) const;
    Param* LocateParam(const string &ParamName);

//fait un vecteur de class Param
    const Param* LocateOrAddParam(const string &ParamName);
    Param* LocateOrAddNonConstParam(const string &ParamName);

//Elimine le param qui � le nom ParamName. Retourne true si trouv� et false sinon
    bool DeleteParamIfExists(const string &ParamName);

// return the number of variable, scanned, parameters (not the fixed ones)
    size_t VariableParamCount() const;


    //Sort dans "stream" le message et les param�tres (cf dump):val, val_t0
    void dump(ostream &stream = cout, const string &Message = "",bool nice_formated_output=false) const;

    //Sort dans "stream" les param�tres: nom val
    void dump_model(ostream &stream = cout, const string &Message = "") const;


    friend ostream& operator << (ostream &stream, const FitParams &p)
    {
        p.dump(stream);
        return stream;
    }

    /*! if FitParam has a name, read the block with the same name else reads the first block and assigns name */
// L'ordre dans le fichier doit �tre name_scanned, minvalue, maxvalue, nombrestep, bool_is_scann, bool_is_time_dependent, tau_var, erro
    void read(const std::string &FileName);
//    void read(istream &S);

// Initialise les param�tres selon la valeur non scann�e de la dataCard (cf datacards.h) ou � la valeur min si scann�e
    void init_value(const std::string &FileName);

// Modification des param�tres
// Nb_steps = 1 signifie que l'on va prendre 2 valeurs min et max. nb_step=2 on prend 3 valeurs: min, moiti� et max  .
// On retourne la condition de fin: false s'il faut encore continuer � modifier et true si on a fini les boucles.
// C'est � dire � priori on ne modifie qu'un param�tre � la fois dans l'ordre d'apparition dans la liste des scans
    bool  Scan_Param();


// Modification des param�tres al�atoirement
    bool  Scan_Param_aleatoire(const gsl_rng * r);

    ~FitParams()
    {
        for (FitParams::iterator i=begin(); i != end(); ++i) delete *i;
    }

    void operator = (const FitParams &);
};

typedef FitParams::iterator ParamIterator;
typedef FitParams::const_iterator ParamCIterator;


vector<string> FitParamNames(const string &FileName);

#endif /* PARAM__H */
