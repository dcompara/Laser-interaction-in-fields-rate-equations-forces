#include <iomanip>
#include <cmath>

#include <cstring> // for strcmp

#include "params.h"


// Sortie du param�tre (val, minv , ...)
void Param::dump(ostream &stream, bool nice_formated_output) const
{
    stream << setprecision(12) << Name() << ' ' << val << ' ' <<  minv << ' ' <<   maxv << ' ' << nbstep << ' ' << boolalpha << is_scanned << ' ' << boolalpha  << is_time_dependent <<  ' ' << tau <<  ' ' << val_t0 ;
}



//pour atteindre un param Retourne le pointeur Param* vers le param qui � le nom ParamName
const Param* FitParams::LocateParam(const string &ParamName) const
{
    ParamCIterator i = begin();
    for ( ; i != end(); ++i)
        if ((*i)->Name() == ParamName) break;
    if (i == end())
    {
        cerr << " Param " << ParamName << " does not exist " << endl;
        return NULL;
    }
    return (*i);
}

//Elimine le param qui � le nom ParamName. Retourne true si trouv� et false sinon
bool FitParams::DeleteParamIfExists(const string &ParamName)
{
    ParamIterator i = begin();
    for ( ; i != end(); ++i)
        if ((*i)->Name() == ParamName)
            {
                erase(i);
                return true;
            }
    return false;
}


//pour atteindre un param Retourne le pointeur Param* vers le param qui � le nom ParamName
Param* FitParams::LocateParam(const string &ParamName)
{
    ParamIterator i = begin();
    for ( ; i != end(); ++i)
        if ((*i)->Name() == ParamName) break;
    if (i == end()) return NULL;
    return (*i);
}

//fait un vecteur const de class Param
const Param * FitParams::LocateOrAddParam(const string &ParamName)
{
    const Param* i = LocateParam(ParamName);
    if (i != NULL) return i;
// no there : add one
    push_back(new Param(ParamName));
    return LocateParam(ParamName);
//  ParamCIterator last = end();
// return (--last); // i.e. the one just inserted (since end() is beyond last)
}

//fait un vecteur de class Param
Param * FitParams::LocateOrAddNonConstParam(const string &ParamName)
{
    Param* i = LocateParam(ParamName);
    if (i != NULL) return i;
// no there : add one
    push_back(new Param(ParamName));
    return LocateParam(ParamName);
//  ParamCIterator last = end();
// return (--last); // i.e. the one just inserted (since end() is beyond last)
}

// return the number of variable (not the fixed ones) parameters
size_t FitParams::VariableParamCount() const
{
    size_t npars = 0;
// count and locate actually variable parameters
    for (ParamCIterator i=begin(); i != end(); ++i)
        if ((*i)->is_scanned) npars++;
    return npars;
}


//Sort dans "stream" le message et les param�tres (cf dump)
void FitParams::dump(ostream &stream, const string &Message, bool nice_formated_output) const
{
    stream << "BEGIN_OF_FITPARAMS " << name << endl;
    if (Message != "") stream << Message << std::endl;

    stream << "Name" << ' ' << "val" << ' ' <<  "minv" << ' ' <<   "maxv" << ' ' << "nbstep" << ' ' << "is_scanned" << ' ' << "is_time_dependent" <<  ' ' << "tau" <<  ' ' << "val_t0" << endl;

    for (ParamCIterator i=begin(); i != end(); ++i)
        {
            (*(*i)).dump(stream,nice_formated_output);
            stream << endl ;
        }

    stream << "END_OF_FITPARAMS " << name << endl;
}

//Sort dans "stream" les param�tres: nom val
void FitParams::dump_model(ostream &stream, const string &Message) const
{
    if (Message != "") stream << Message << std::endl;
    for(size_t i=0; i< this->size(); ++i)
        {
            stream  << setprecision(12) << (*this)[i]->name << " " <<(*this)[i]->val << endl  ;
        }
}




#include <fstream>

#include "datacards.h"
//permet de recuperer les diff�rents Param dans un vector de Param
// La ligne contient
// Le nom, Valeur_in, Valeur_fin, Nb_steps,  is_scanned ("true" or "false"),  is_time_dependent("true" or "false"), tau;
void FitParams::read(const std::string &FileName)
{
    DataCards data_list(FileName.c_str()); // Lit le fichier et cr�e les datacards.
    double Tau_Modif = data_list.DParam("Tau_Modif"); // Valuer par d�faut de la variation temporelle

    ifstream S(FileName.c_str());

    string line;
    while (getline(S, line))
        {
            if(line == "BEGIN_OF_FITPARAMS") break;
        }

    while (getline(S, line))
        {
            if (line == "END_OF_FITPARAMS") break;

            string name, boolean; // nom de la variable et nom avec SCAN_
            if (line.compare(0,6,"@SCAN_") == 0)  // C'est un param�tre a scanner
                {
                    size_t p = 6; // longueur de "@SCAN_"
                    size_t q = line.find_first_of(" \t"); // Prend la fin du nom (s�par� des autres par TAB)
                    name = line.substr(p,q-p); // on enl�ve @SCAN_ pour prendre le nom
                }

            Param *p=LocateOrAddNonConstParam(name); // AJOUTE UN PARAMETRE DU NOM Name
            char name_scanned [100];
            char bool_is_scann [100]= {'f','a','l','s','e'};
            char bool_is_time_dependent [100]= {'f','a','l','s','e'};
            p->tau  = Tau_Modif; //Met par d�faut la variation temporelle de Tau_Modif. Sera modifier ci dessous si une autre valeur est donn�e

// Lit le fichier. Si il n'y a pas toute les donn�es il n'affecte pas les valeurs sui restent � leur valeur par d�faut
            std::sscanf (line.c_str(),"%s %lf %lf %d %s %s %lf %lf", name_scanned, &(p->minv),  &(p->maxv), &(p->nbstep), bool_is_scann, bool_is_time_dependent, &(p->tau), &(p->val_t0));
            boolean = bool_is_scann;
            if (boolean == "true")  p->is_scanned = true;
            else p->is_scanned = false; // is_scanned(false) par d�faut
            boolean = bool_is_time_dependent;
            if (boolean == "true")  p->is_time_dependent = true;
            else  p->is_time_dependent = false; // is_time_dependent(false) par d�faut
        }

    S.close();
}

//Lit les param�tres dans le "stream"
//void FitParams::read(istream &S)
//{
//
//}


void FitParams::operator = (const FitParams &other)
{


//  if (name == "NoName" ) name = other.name;
// don't know why the name was not copied.
// A real copy constructor should copy the name as well:
    name = other.name;

// first clear the parameters in this instance that are not in the other
    FitParams::iterator i=begin();
    while(i != end())
        {
            if ( ! other.LocateParam((*i)->name))
                {
                    delete *i; // free memory (delete parameters)
                    i = erase(i);
                }
            else
                {
                    ++i;
                }
        }


    for (FitParams::const_iterator i=other.begin(); i != other.end(); ++i)
        {
            Param *p = LocateOrAddNonConstParam((*i)->name); // get a pointer to the parameter of this instance
            *p = *(*i); // copy the content of parameters
        }

}

FitParams::FitParams(const FitParams &other) :
    vector<Param*>()
{
    (*this) = other;
}

#include <cstdio>

vector<string> FitParamNames(const string &FileName)
{
    vector<string> names;
    FILE* f = fopen(FileName.c_str(), "r");

    char line[4096];
    while (fgets(line,4096,f))
        {
            char name[100];
            if (sscanf(line,"BEGIN_OF_FITPARAMS %s",name) == 1)
                names.push_back(name);
        }
    return names;
}

// Initialise les param�tres selon la valeur non scann�e de la dataCard (cf datacards.h) ou � la valeur min si scann�e
void FitParams::init_value(const std::string &FileName)
{
    DataCards data_list(FileName.c_str()); // Lit le fichier et cr�e les datacards.
    for (DataCards::CardList::iterator i = data_list.cards.begin(); i != data_list.cards.end(); i++)
        Param *p=LocateOrAddNonConstParam((*i).kw ); // AJOUTE (si n'existe pas) UN PARAMETRE DU NOM du param�tre dans la dataCard

    for (ParamIterator i=(*this).begin(); i != (*this).end(); ++i) // boucle sur les param�tres
        {
            Param &p = **i;

            if( p.is_scanned == false && data_list.HasKey(p.name)) // Si param�tre ne sera pas � scanner
                {
                    p.val = data_list.DParam(p.name); // Si param�tre constant on prend sa valeur dans la datacard(si elle existe)
                    p.val_t0 = p.val;
                }
            else
                {
                    p.val =  p.minv;           // Si param�tre � scanner on l'initialise � sa valeur de d�part (min)
                    p.val_t0 = p.val;
                }

        }

}




// Modification des param�tres dans l'ordre
// Nb_steps = 1 signifie que l'on va prendre 2 valeurs min et max. nb_step=2 on prend 3 valeurs: min, moiti� et max  .
// On retourne la condition de fin: false s'il faut encore continuer � modifier et true si on a fini les boucles.
// C'est � dire � priori on ne modifie qu'un param�tre � la fois dans l'ordre d'apparition dans la liste des scans
bool  FitParams::Scan_Param()
{
    for (ParamIterator i=(*this).begin(); i != (*this).end(); ++i) // boucle sur les param�tres
        {
            Param &p = **i;
            // cout << p << endl << endl;
            if( p.is_scanned == false ) // Pas de boucle pour ce param�tre
                continue;

            if(p.val_t0 >= p.maxv) // On a fini la i�me boucle on r�initailise ce param�tre et on passe � la boucle suivante
                {
                    p.val_t0  = p.minv;
                    p.val = p.val_t0;
                }
            else // On n'atteind pas encore la fin de la boucle
                {
                    p.val_t0 += (p.maxv - p.minv)/p.nbstep; // On augmente le param�tre
                    p.val = p.val_t0;
                    return false;                       // et on arrete donc les boucles
                }
        }
    return true;
}


// Modification des param�tres al�atoirement
bool  FitParams::Scan_Param_aleatoire(const gsl_rng * r)
{
    for (ParamIterator i=(*this).begin(); i != (*this).end(); ++i) // boucle sur les param�tres
        {
            Param &p = **i;

            if( p.is_scanned == false ) // Pas de boucle pour ce param�tre
                continue;
            else
                {
                    p.val_t0 = gsl_ran_flat (r, p.minv, p.maxv); // tirage au hasard entre les valeur min et max possibles
                    p.val = p.val_t0;
                }
            cout << p.name << " " << p.val_t0 <<  endl;
        }
    return false;
}



