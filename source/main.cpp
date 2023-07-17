#include "ReadAnyBURLRules.cpp"

using namespace std;

int main(int argc, char* argv[]){
    std::string dataset(argv[1]);
    std::string alg(argv[2]);
    std::string topr(argv[3]);
    std::string beta(argv[4]);
    std::string len(argv[5]);

    if(alg == "ExactD+")  // finished
        Effective_hg_global_f1(dataset, std::stod(topr), "df1");
    if(alg == "ExactH+")  // finished
        Effective_hg_global_f1(dataset, std::stod(topr), "hf1");
    else if(alg == "ExactD")  // finished
        Effective_hg_global_greater_f1(dataset, std::stod(topr), "df1");
    else if(alg == "ExactH")  // finished
        Effective_hg_global_greater_f1(dataset, std::stod(topr), "hf1");
    else if(alg == "GloD")  // finished
        Effective_prop_opt_global_cross(dataset, topr, "df1");
    else if(alg == "GloH")  // finished
        Effective_prop_opt_global_cross(dataset, topr, "hf1");
    else if(alg == "PerD")  // finished
        Effective_prop_opt_personalized_cross(dataset, topr, std::stod(beta), "d");
    else if(alg == "PerH")  // finished
        Effective_prop_opt_personalized_cross(dataset, topr, std::stod(beta), "h");
    else if(alg == "PerD+")  // finished
        Effective_prop_opt_personalized_cross(dataset, topr, std::stod(beta), "dp");
    else if(alg == "PerH+")  // finished
        Effective_prop_opt_personalized_cross(dataset, topr, std::stod(beta), "hp");

    //else if(alg == "mpcount")
    //    MetaPathCount(dataset, 1);
    //else if(alg == "lensplit")
    //    MetaPathLenSplit(dataset);
    //else if(alg == "hg_stats")
    //    Effective_hg_stats(dataset);
    //else if(alg == "global_scale")
    //    Scalability_prop_opt_global(dataset, topr, method, len);
    //else if(alg == "personalized_scale")
    //    Scalability_prop_opt_personalized(dataset, topr, std::stod(beta), method, len);
    //else if(alg == "hg_scale")
    //    Scalability_hg_greater_f1(dataset, std::stod(topr), method, len);
    //else if(alg == "matching_graph_time")
    //    MatchingGraphTime(dataset);
    //else if(alg == "union")
    //    Effective_hg_global_f1_by_union(dataset, std::stod(topr), method);
}