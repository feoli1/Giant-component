#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <array>
#include <vector>
#include <memory>
#include <iomanip>
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>
#include <stdexcept>

std::random_device rd;
std::mt19937 generateur(rd());

// Problèmes à régler
// Taille des composantes, ligne ou totale
// Mauvaise utilisation de la cache pour générer une ligne (>1T)

class Reseau
{
    struct Node
    {
        double x;
        double y;
        Node* parent;
        uint_least32_t lignesRestantes;
        uint_least64_t taille;
        Node* origine;
    };
    struct Ligne
    {
        std::vector<Node> nodes;
        std::unique_ptr<uint_least32_t[]> chunkId;
        std::vector<Node*> touche;
    };
public:
    Reseau(uint_least64_t N, double densite, double portee=1.)
    {
        portee2 = portee*portee;
        double tailleMonde = std::sqrt(N/densite);
        nbChunks = tailleMonde / portee; // entier, arrondi vers le bas -> les chunks sont au moins la longeur de la portee
        if (nbChunks < 2)
        {
            std::cout << "Le monde est trop petit!";
            std::exit(1); // Meilleure solution?
        }
        tailleChunk = tailleMonde / nbChunks; //double
        noeudsRestants = N;
        lignesRestantes = nbChunks;
        randDouble = std::uniform_real_distribution<double>(0., 1.);
        randChunk = std::uniform_int_distribution<uint_least32_t>(1, nbChunks);
        baseLigne.chunkId = std::unique_ptr<uint_least32_t[]>(new uint_least32_t[nbChunks+1]);
        ancLigne.chunkId  = std::unique_ptr<uint_least32_t[]>(new uint_least32_t[nbChunks+1]);
        novLigne.chunkId  = std::unique_ptr<uint_least32_t[]>(new uint_least32_t[nbChunks+1]);
        baseLigne.chunkId[0] = 0;
        ancLigne.chunkId[0] = 0;
        novLigne.chunkId[0] = 0;
        genererLigne(baseLigne);
        genererLigne(novLigne);
        maxTaille = 0;
        générerComposantes(baseLigne, novLigne);
        for (Node* origine = baseLigne.nodes.data(); origine < baseLigne.nodes.data() + baseLigne.chunkId[nbChunks]; origine++)
        {
            origine->parent = racine(origine);
        }
        for (Node* origine = baseLigne.nodes.data(); origine < baseLigne.nodes.data() + baseLigne.chunkId[nbChunks]; origine++)
        {
            Node* rac = racine(origine);
            if (rac->lignesRestantes == lignesRestantes)
            {
                if (!rac->origine) // La racine est touchée pour la première fois, on lui donne une origine
                {
                    rac->origine = origine;
                    novLigne.touche.push_back(rac);
                }
                else origine->parent = rac->origine; // Combine les composantes par la base
            }
        }
        while (lignesRestantes > 0)
        {
            std::swap(novLigne, ancLigne);
            genererLigne(novLigne);
            générerComposantes(ancLigne, novLigne);
            // Gestion de la ligne de base
            uint_least32_t iMax = ancLigne.touche.size(); //On initialise la taille au préalable car la taille peut augementer
            for (uint_least32_t i  = 0; i < iMax; i++)
            {
                Node* nodeTouche = ancLigne.touche[i];
                Node* origine = nodeTouche->origine;
                Node* rac = racine(nodeTouche);
                if (rac->lignesRestantes > lignesRestantes) // La racine est restée sur l'ancienne ligne, fin de la composante
                {
                    if (!rac->origine) // La racine est touchée pour la première fois, on lui donne une origine
                    {
                        rac->origine = origine;
                        ancLigne.touche.push_back(rac); // Pour le nettoyage de l'ancienne ligne
                    }
                    origine->taille = rac->taille; // NE PAS METTRE DANS LE IF, efficace même si la racine n'a pas bougé lorsque la nouvelle composante a été générée. On pourrait remplacer origine->taille par rac->origine->taille mais l'effet est le même
                    origine->parent = rac->origine; // Retourne la racine à l'origine OU combine les composantes par la base
                }
                else
                {
                    if (!rac->origine) // La racine est touchée pour la première fois, on lui donne une origine
                    {
                        rac->origine = origine;
                        novLigne.touche.push_back(rac);
                    }
                    else origine->parent = rac->origine; // Combine les composantes par la base
                }
            }
            for (Node* nodeTouche : ancLigne.touche)
            {
                nodeTouche->origine = NULL; // Réinitialisation de l'ancienne ligne
            }
            ancLigne.touche.clear(); // Réinitialisation de l'ancienne ligne
        }
        for (Node* nodeTouche : novLigne.touche)
        {
            nodeTouche->origine->parent = nodeTouche; // Les parents de la ligne de base sont mis à jour
        }
        générerComposantes(novLigne, baseLigne);
        GCFraction = double(maxTaille) / N;
    }
    double GCFraction;

private:
    Node* racine(Node* i) // Un noeud a toujours
    {
        while(i != i->parent)
        {
            i->parent = i->parent->parent;
            i = i->parent;
        }
        return i;
    }
    
    void unir(Node* p, Node* q)
    {
        Node* i = racine(p);
        Node* j = racine(q);
        if (i != j) // Les deux ne sont pas déjà dans la même composante
        {
                 if (i->lignesRestantes > j->lignesRestantes) { i->parent = j; j->taille += i->taille; if (j->taille > maxTaille) maxTaille = j->taille; } // La racine doit toujours être transférée à la nouvelle ligne
            else if (j->lignesRestantes > i->lignesRestantes) { j->parent = i; i->taille += j->taille; if (i->taille > maxTaille) maxTaille = i->taille; }
            else if (i->taille > j->taille)                   { j->parent = i; i->taille += j->taille; if (i->taille > maxTaille) maxTaille = i->taille; } // Si les deux racines sont sur la même ligne, la composante la plus grande l'emporte
            else                                              { i->parent = j; j->taille += i->taille; if (j->taille > maxTaille) maxTaille = j->taille; }
        }
    }

    void genererLigne(Ligne& ligne)
    {
        uint_least32_t noeudsLigne;
        if (lignesRestantes == 1) noeudsLigne = noeudsRestants;
        else
        {
            double p = 1.0 / lignesRestantes;
            double q = 1 - p;
            noeudsLigne = noeudsRestants*p + sqrt(2*noeudsRestants*p*q) * boost::math::erf_inv(2*randDouble(generateur)-1);
            if (noeudsLigne > noeudsRestants) noeudsLigne = noeudsRestants; // Statistiquement impossible pour de grands mondes
        }
        lignesRestantes--;
        noeudsRestants -= noeudsLigne;

        uint_least32_t* chunkId = ligne.chunkId.get();
        // On calcule le pointeur vers le début de chaque chunk
        std::fill(chunkId+1, chunkId+nbChunks+1, 0); // Réinitialisation
        for (uint_least32_t i = 0; i < noeudsLigne; i++) // On distribue les noeuds aux chunks
        {
            chunkId[randChunk(generateur)]++;
        }
        for (uint_least32_t i = 1; i < nbChunks; i++) // Le premier noeud du chunk i est chunkId[i], le passé le dernier est chunkId[i+1]
        {
            chunkId[i+1] += chunkId[i];
        }
        // On calcule les noeuds
        if (noeudsLigne > ligne.nodes.capacity()) // On s'assure que l'array est suffisemment grand. On évite également que les anciennes données soient déplacées
        {
            ligne.nodes.clear(); // Clear pour ne pas déplacer la mémoire
            uint_least32_t noeudsLigneMarge = 1.1*noeudsLigne; // Marge de manoeuvre, un peu plus de mémoire utilisée mais beaucoup moins de réallocations
            ligne.nodes.reserve(noeudsLigneMarge);
            for (Node* node = ligne.nodes.data(); node < ligne.nodes.data() + noeudsLigneMarge; node++)
            {
                node->origine = NULL;
            }
        }
        for (uint_least32_t i = 0; i < noeudsLigne; i++)
        {
            ligne.nodes[i] = { tailleChunk * randDouble(generateur),
                               tailleChunk * randDouble(generateur),
                               &ligne.nodes[i],
                               lignesRestantes,
                               1 };
        }
    }
    
    void générerComposantes(Ligne& ancLigne, Ligne& novLigne)
    {
        for (uint_least32_t k = 0; k < nbChunks; k++) // On itère sur tous les chunks
        {
            for (Node* p = ancLigne.nodes.data() + ancLigne.chunkId[k]; p < ancLigne.nodes.data() + ancLigne.chunkId[k+1]; p++)
            {
                // Chunk central
                for (Node* q = p + 1; q < ancLigne.nodes.data() + ancLigne.chunkId[k+1]; q++)
                {
                    double dx = q->x - p->x;
                    double dy = q->y - p->y;
                    if (dx*dx + dy*dy < portee2) unir(p, q);
                }
                // Chunk est
                uint_least32_t l;
                if (k == nbChunks - 1) l = 0; else l = k + 1;
                for (Node* q = ancLigne.nodes.data() + ancLigne.chunkId[l]; q < ancLigne.nodes.data() + ancLigne.chunkId[l+1]; q++)
                {
                    double dx = q->x + tailleChunk - p->x;
                    double dy = q->y - p->y;
                    if (dx*dx + dy*dy < portee2) unir(p, q);
                }
                // Chunk sud
                for (Node* q = novLigne.nodes.data() + novLigne.chunkId[k]; q < novLigne.nodes.data() + novLigne.chunkId[k+1]; q++)
                {
                    double dx = q->x - p->x;
                    double dy = q->y + tailleChunk - p->y;
                    if (dx*dx + dy*dy < portee2) unir(p, q);
                }
                // Chunk sud-ouest
                if (k == 0) l = nbChunks - 1; else l = k - 1;
                for (Node* q = novLigne.nodes.data() + novLigne.chunkId[l]; q < novLigne.nodes.data() + novLigne.chunkId[l+1]; q++)
                {
                    double dx = q->x - tailleChunk - p->x;
                    double dy = q->y + tailleChunk - p->y;
                    if (dx*dx + dy*dy < portee2) unir(p, q);
                }
                // Chunk sud-est
                if (k == nbChunks - 1) l = 0; else l = k + 1;
                for (Node* q = novLigne.nodes.data() + novLigne.chunkId[l]; q < novLigne.nodes.data() + novLigne.chunkId[l+1]; q++)
                {
                    double dx = q->x + tailleChunk - p->x;
                    double dy = q->y + tailleChunk - p->y;
                    if (dx*dx + dy*dy < portee2) unir(p, q);
                }
            }
        }
    }

    double portee2;
    uint_least32_t nbChunks;
    std::uniform_real_distribution<double> randDouble;
    std::uniform_int_distribution<uint_least32_t> randChunk;
    double tailleChunk;
    Ligne baseLigne;
    Ligne ancLigne;
    Ligne novLigne;
    uint_least64_t noeudsRestants;
    uint_least32_t lignesRestantes;
    uint_least64_t maxTaille;
};

int main(int argc, char *argv[])
{
    uint_least64_t N;
    double densite;
    int simulations;
    std::unique_ptr<std::string[]> Argv(new std::string[argc]);
    for (int i = 0; i < argc; i++)
    {
        Argv[i] = std::string(argv[i]);
    }
    std::string* p;
    p = std::find(Argv.get(), Argv.get() + argc, "-n");
    if (p == Argv.get() + argc | p == Argv.get() + argc - 1) throw std::logic_error("Aucun nombre de noeuds fourni (-n)");
    N = std::stoull(*(p+1));
    p = std::find(Argv.get(), Argv.get() + argc, "-d");
    if (p == Argv.get() + argc | p == Argv.get() + argc - 1) throw std::logic_error("Aucune densité fournie (-d)");
    densite = std::stod(*(p+1));
    p = std::find(Argv.get(), Argv.get() + argc, "-i");
    if (p == Argv.get() + argc | p == Argv.get() + argc - 1) throw std::logic_error("Aucun nombre de simulations fourni (-i)");
    simulations = std::stoi(*(p+1));
    
    std::cout << N << ", densité: " << densite << ", simulations: " << simulations << std::endl;
    std::unique_ptr<double[]> fractions(new double[simulations]);
    std::ostringstream folder;
    folder << "Data/" << densite;
    std::filesystem::create_directories(folder.str());
    std::ofstream file;
    folder << "/" << N;
    file.open(folder.str(), std::ios_base::app);
    for (int i = 0; i < simulations; i++)
    {
        std::cout << i << "/" << simulations << " simulations compétées...\n";
        double fraction = Reseau(N, densite).GCFraction;
        std::cout << fraction << std::endl;
        file << N << " " << densite << " " << fraction << std::endl;
        fractions[i] = fraction;
    }
    std::cout << simulations << "/" << simulations << " simulations compétées.\n";
    double mean = 0;
    for (int i = 0; i < simulations; i++)
    {
        mean += fractions[i];
    }
    mean /= simulations;
    double err = 0;
    for (int i = 0; i < simulations; i++)
    {
        err += (fractions[i]-mean)*(fractions[i]-mean);
    }
    err = std::sqrt(err) / simulations; // un écart type
    std::cout << "Fraction:" << mean << " ± " << err << std::endl;
    return 0;
}
