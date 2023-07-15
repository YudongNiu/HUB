#include "ReadAnyBURLRules.cpp"

using namespace std;

void dblp_coauthor(){
    auto qp = new Pattern();
    qp->ETypes.push_back(0);
    qp->NTypes.push_back(-1);
    qp->NTypes.push_back(-1);
    qp->EDirect.push_back(1);

    const unsigned int asso = 1;

    HeterGraph g("DBLP");

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();

    while(visited->size()< qp->ETypes.size()+1){
        visited->push_back(new std::vector<bool>(g.NT.size(), false));
        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
    }

    auto frontiers = new std::vector<unsigned int>();
    auto peers = new std::vector<unsigned int>();
    Peers(qp, &g, frontiers, peers, visited);

    std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, &g, visited, back_visited);
    auto ractive = new std::vector<std::vector<unsigned int> *>();
    for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g.NT.size(), 0));
    for (unsigned int n = 0; n < g.NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

    auto hidden_graph = hidden_graph_construction(asso, qp, &g, visited, ractive);

    ofstream outfile;
    outfile.open("file.dat", ios::out);

    for(unsigned int i=0;i<hidden_graph->size();i++){
        for(unsigned int j=0;j<hidden_graph->at(i)->size();j++){
            if(hidden_graph->at(i)->at(j) != i)
                outfile<<i<<" "<<hidden_graph->at(i)->at(j)<<std::endl;
        }
    }
    outfile.close();
}

int main(int argc, char* argv[]){
    /** cod */
    std::string dataset(argv[1]);
    std::string method(argv[2]);
    std::string topr(argv[3]);
    std::string beta(argv[4]);
    std::string len(argv[5]);
    std::string experiment(argv[6]);

    if(experiment == "effect_hg_global")  // finished
        Effective_hg_global_f1(dataset, std::stod(topr), method);
    else if(experiment == "effect_hg_global_greater")  // finished
        Effective_hg_global_greater_f1(dataset, std::stod(topr), method);
    else if(experiment == "effect_prop_personalized_cross")  // finished
        Effective_prop_opt_personalized_cross(dataset, topr, std::stod(beta), method);
    else if(experiment == "effect_prop_global_cross")  // finished
        Effective_prop_opt_global_cross(dataset, topr, method);
    else if(experiment == "mpcount")
        MetaPathCount(dataset, 1);
    else if(experiment == "lensplit")
        MetaPathLenSplit(dataset);
    else if(experiment == "hg_stats")
        Effective_hg_stats(dataset);
    else if(experiment == "global_scale")
        Scalability_prop_opt_global(dataset, topr, method, len);
    else if(experiment == "personalized_scale")
        Scalability_prop_opt_personalized(dataset, topr, std::stod(beta), method, len);
    else if(experiment == "hg_scale")
        Scalability_hg_greater_f1(dataset, std::stod(topr), method, len);
    else if(experiment == "matching_graph_time")
        MatchingGraphTime(dataset);
    else if(experiment == "union")
        Effective_hg_global_f1_by_union(dataset, std::stod(topr), method);

    /** density filter */
//    std::string dataset(argv[1]);
//    std::string dens_threshod(argv[2]);
//    ReadAnyBURLRules_DensityFilter(dataset, std::stod(dens_threshod));

    /** dblp-coauthor-extraction */
//    dblp_coauthor();
}