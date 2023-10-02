#include "cod.cpp"
#include "effectiveness.cpp"
#include <chrono>
#include <thread>

std::vector<std::string> split(const std::string& str, const std::string &pattern, bool app){
    std::vector<std::string> res;
    if(str == "") return res;

    std::string::size_type temp_pos;
    std::string appstr;
    if(app) appstr = str + pattern;
    else appstr = str;

    int size = appstr.size();
    for(unsigned int i=0;i<size;i++){
        temp_pos = appstr.find(pattern, i);
        if(temp_pos < size){
            std::string substr = appstr.substr(i, temp_pos - i);
            res.push_back(substr);
            i = temp_pos;
        }
    }
    return res;
}

//void ReadAnyBURLRules_DensityFilter(const std::string &choice, double dens_threshold){
//    HeterGraph g(choice);
//
//    std::string rules_path = (choice + "/cod-gnn-rules.dat");
//    std::ifstream rules_in;
//    rules_in.open(rules_path);
//    std::string rules_line;
//
//    unsigned int rule_count = 0;
//
//    auto visited = new std::vector<std::vector<bool>*>();
//    auto back_visited = new std::vector<std::vector<bool>*>();
//    auto qp = new Pattern();
//
//    getline(rules_in, rules_line);
//    int state = -1;
//    // state == 0: next int represents a variable rule;
//    // state == 1: next int represents a instance rule;
//
//    std::string::size_type temp_pos;
//    int size = rules_line.size();
//    for(unsigned int i = 0; i<size;i++) {
//        temp_pos = rules_line.find(' ', i);
//        if (temp_pos < size) {
//            int sub = stoi(rules_line.substr(i, temp_pos - i));
//
//            if (sub == -1) state = 0; // next int indicates a variable rule
//            else if (sub == -2) qp->EDirect.push_back(1); // ->
//            else if (sub == -3) qp->EDirect.push_back(-1); // <-
//            else if (sub == -4) { // pop
//                qp->EDirect.pop_back();
//                qp->ETypes.pop_back();
//                qp->NTypes.pop_back();}
//            else if (sub == -5) state = 1; // next int indicates a instance rule
//            else {
//                if (state == 1) { // instance rule
//                    qp->instance = sub;
//                    while(visited->size()< qp->ETypes.size()+1){
//                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
//                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
//                    }
//                    prop::density(qp, &g, visited, back_visited, dens_threshold);
//                    qp->instance = -1;
//                    state = -1;
//                } else if (state == 0) { // variable rule
//                    qp->ETypes.push_back(sub);
//                    qp->NTypes.push_back(-1);
//                    while(visited->size()< qp->ETypes.size()+1){
//                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
//                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
//                    }
//                    prop::density(qp, &g, visited, back_visited, dens_threshold);
//                    state = -1;
//                } else {
//                    qp->ETypes.push_back(sub);
//                    qp->NTypes.push_back(-1);
//                }
//            }
//            i = temp_pos;
//        }
//    }
//    rules_in.close();
//
//    for (auto &it : *visited) delete it;
//    for (auto &it : *back_visited) delete it;
//    delete visited;
//    delete back_visited;
//}

void MetaPathLenSplit(const std::string &choice){
    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    std::string rules_path_len1 = (choice + "/"+choice+"-cod-global-rules.dat1");
    std::ofstream rules_out1;
    rules_out1.open(rules_path_len1);

    std::string rules_path_len2 = (choice + "/"+choice+"-cod-global-rules.dat2");
    std::ofstream rules_out2;
    rules_out2.open(rules_path_len2);

    std::string rules_path_len3 = (choice + "/"+choice+"-cod-global-rules.dat3");
    std::ofstream rules_out3;
    rules_out3.open(rules_path_len3);

    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    int state = -1;
    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;

                    unsigned int len = qp->ETypes.size() - 1;
                    if(len == 1){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(qp->EDirect[jj] == 1) rules_out1 << "-2 ";
                            else rules_out1<< "-3 ";
                            rules_out1<<qp->ETypes[jj]<<" ";
                        }
                        rules_out1<<"-5 "<<qp->instance<<" ";
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out1<<"-4 ";
                    }
                    else if(len == 2){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(qp->EDirect[jj] == 1) rules_out2 << "-2 ";
                            else rules_out2<< "-3 ";
                            rules_out2<<qp->ETypes[jj]<<" ";
                        }
                        rules_out2<<"-5 "<<qp->instance<<" ";
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out2<<"-4 ";
                    }
                    else if(len == 3){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(qp->EDirect[jj] == 1) rules_out3 << "-2 ";
                            else rules_out3<< "-3 ";
                            rules_out3<<qp->ETypes[jj]<<" ";
                        }
                        rules_out3<<"-5 "<<qp->instance<<" ";
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out3<<"-4 ";
                    }

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);

                    unsigned int len = qp->ETypes.size();
                    if(len == 1){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(jj == qp->ETypes.size()-1) rules_out1 <<"-1 ";
                            if(qp->EDirect[jj] == 1) rules_out1 << "-2 ";
                            else rules_out1<< "-3 ";
                            rules_out1<<qp->ETypes[jj]<<" ";
                        }
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out1<<"-4 ";
                    }
                    else if(len == 2){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(jj == qp->ETypes.size()-1) rules_out2 <<"-1 ";
                            if(qp->EDirect[jj] == 1) rules_out2 << "-2 ";
                            else rules_out2<< "-3 ";
                            rules_out2<<qp->ETypes[jj]<<" ";
                        }
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out2<<"-4 ";
                    }
                    else if(len == 3){
                        for(int jj=0;jj<qp->ETypes.size();jj++){
                            if(jj == qp->ETypes.size()-1) rules_out3 <<"-1 ";
                            if(qp->EDirect[jj] == 1) rules_out3 << "-2 ";
                            else rules_out3<< "-3 ";
                            rules_out3<<qp->ETypes[jj]<<" ";
                        }
                        for(int jj=0;jj<qp->ETypes.size();jj++) rules_out3<<"-4 ";
                    }
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    rules_out1.close();
    rules_out2.close();
    rules_out3.close();
}

void MetaPathCount(const std::string &choice, unsigned int limited){
    unsigned int LEN = 5;

    std::string rules_path;
    if (limited == 0) rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    else rules_path = (choice + "/"+choice+"-cod-global-rules.limit");

    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto qp = new Pattern();
    unsigned int rule_count = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    std::vector<unsigned int> path_len_count;
    for(unsigned int i=0;i<=LEN;i++) path_len_count.push_back(0);

    int state = -1;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                    qp->EDirect.pop_back();
                    qp->ETypes.pop_back();
                    qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;

                    unsigned int len = qp->ETypes.size() - 1;
                    if(len < LEN) path_len_count[len] += 1;

                    rule_count++;
                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);

                    unsigned int len = qp->ETypes.size();
                    if(len < LEN) path_len_count[len] += 1;

                    rule_count++;
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    std::cout<<"Length Count:"<<std::endl;
    for(unsigned int i=1;i<=LEN;i++)
        std::cout<<"len="<<i<<"  count:"<<path_len_count[i]<<std::endl;
    std::cout<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void MatchingGraphTime(const std::string &choice){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    bool found = effectiveness::COD_matching_graph_time(qp, &g, visited, back_visited, avg_time);
                    if(found) rule_count++;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool found = effectiveness::COD_matching_graph_time(qp, &g, visited, back_visited, avg_time);
                    if(found) rule_count++;
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000;

    std::cout<<"~matching_graph_time_per_rule:"<<avg_time<<" ms"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Scalability_prop_opt_global(const std::string &choice, const std::string &topr, const std::string &method, const std::string &len){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat"+len);
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    bool found = effectiveness::COD_prop_global_scale(qp, &g, std::stod(topr), visited, back_visited, method, avg_time);
                    if(found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool found =  effectiveness::COD_prop_global_scale(qp, &g, std::stod(topr), visited, back_visited, method, avg_time);
                    if(found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_prop_opt_global_cross(const std::string &choice, const std::string &topr, const std::string &method){
    std::ifstream hg_dom_greater_in, hg_dom_in;
    if(method == "df1") {
        hg_dom_greater_in.open("global_res/" + choice + "/df1/hg_global_greater_r"+topr+".res");
        hg_dom_in.open("global_res/"+choice+"/df1/hg_global_r"+topr+".res");
    }
    else {
        hg_dom_greater_in.open("global_res/" + choice + "/hf1/hg_global_greater_r"+topr+".res");
        hg_dom_in.open("global_res/"+choice+"/hf1/hg_global_r"+topr+".res");
    }
    std::string hg_dom_line;

    auto hg_doms_greater = new std::vector<std::vector<unsigned int>*>();
    auto hg_isdoms_greater = new std::vector<bool>();

    while(getline(hg_dom_greater_in, hg_dom_line)){
        if(hg_dom_line[0] == ')' || hg_dom_line[0] == '('){
            auto temp_doms = split(hg_dom_line, " ", false);
            hg_isdoms_greater->push_back(hg_dom_line[0] == '(');
            hg_doms_greater->push_back(new std::vector<unsigned int>());
            for(unsigned int i=1;i<temp_doms.size();i++){
                hg_doms_greater->at(hg_doms_greater->size()-1)->push_back((unsigned int)stoi(temp_doms[i]));
            }
        }
    }
    hg_dom_greater_in.close();

    auto hg_doms = new std::vector<std::vector<unsigned int>*>();
    auto hg_isdoms = new std::vector<bool>();

    while(getline(hg_dom_in, hg_dom_line)){
        if(hg_dom_line[0] == ')' || hg_dom_line[0] == '('){
            auto temp_doms = split(hg_dom_line, " ", false);
            hg_isdoms->push_back(hg_dom_line[0] == '(');
            hg_doms->push_back(new std::vector<unsigned int>());
            for(unsigned int i=1;i<temp_doms.size();i++){
                hg_doms->at(hg_doms->size()-1)->push_back((unsigned int)stoi(temp_doms[i]));
            }
        }
    }
    hg_dom_in.close();

    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;
    unsigned int dom_count = 0;
    double avg_relative_error = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                    qp->EDirect.pop_back();
                    qp->ETypes.pop_back();
                    qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    double relative_error = effectiveness::COD_prop_global_cross_f1(qp, &g, std::stod(topr), visited, back_visited, method, hg_doms, hg_isdoms, hg_doms_greater, hg_isdoms_greater, dom_count, avg_time);
                    if(relative_error >= -0.1){
                        rule_count++;
                        avg_relative_error += relative_error;
                    }

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    double relative_error = effectiveness::COD_prop_global_cross_f1(qp, &g, std::stod(topr), visited, back_visited, method, hg_doms, hg_isdoms, hg_doms_greater, hg_isdoms_greater, dom_count, avg_time);
                    if(relative_error > -0.1){
                        rule_count++;
                        avg_relative_error += relative_error;
                    }
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_relative_error /= rule_count;
    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~goodness:"<<avg_relative_error<<std::endl;
    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Scalability_prop_opt_personalized(const std::string &choice, const std::string &topr, double beta, const std::string &method, const std::string &len){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat"+len);
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    bool found = effectiveness::COD_prop_personalized_scale(qp, &g, std::stod(topr), beta, visited, back_visited, method, avg_time);
                    if(found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool found = effectiveness::COD_prop_personalized_scale(qp, &g, std::stod(topr), beta, visited, back_visited, method, avg_time);
                    if(found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;

                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_prop_opt_personalized_cross(const std::string &choice, const std::string &topr, double beta, const std::string &method){
    std::ifstream hg_dom_greater_in, hg_dom_in;
    if(method == "df1" || method == "dp") {
        hg_dom_greater_in.open("global_res/" + choice + "/df1/hg_global_greater_r"+topr+".res");
        hg_dom_in.open("global_res/"+choice+"/df1/hg_global_r"+topr+".res");
    }
    else {
        hg_dom_greater_in.open("global_res/" + choice + "/hf1/hg_global_greater_r"+topr+".res");
        hg_dom_in.open("global_res/"+choice+"/hf1/hg_global_r"+topr+".res");
    }
    std::string hg_dom_line;

    auto hg_doms_greater = new std::vector<std::vector<unsigned int>*>();
    auto hg_isdoms_greater = new std::vector<bool>();

    while(getline(hg_dom_greater_in, hg_dom_line)){
        if(hg_dom_line[0] == ')' || hg_dom_line[0] == '('){
            auto temp_doms = split(hg_dom_line, " ", false);
            hg_isdoms_greater->push_back(hg_dom_line[0] == '(');
            hg_doms_greater->push_back(new std::vector<unsigned int>());
            for(unsigned int i=1;i<temp_doms.size();i++){
                hg_doms_greater->at(hg_doms_greater->size()-1)->push_back((unsigned int)stoi(temp_doms[i]));
            }
        }
    }
    hg_dom_greater_in.close();

    auto hg_doms = new std::vector<std::vector<unsigned int>*>();
    auto hg_isdoms = new std::vector<bool>();

    while(getline(hg_dom_in, hg_dom_line)){
        if(hg_dom_line[0] == ')' || hg_dom_line[0] == '('){
            auto temp_doms = split(hg_dom_line, " ", false);
            hg_isdoms->push_back(hg_dom_line[0] == '(');
            hg_doms->push_back(new std::vector<unsigned int>());
            for(unsigned int i=1;i<temp_doms.size();i++){
                hg_doms->at(hg_doms->size()-1)->push_back((unsigned int)stoi(temp_doms[i]));
            }
        }
    }
    hg_dom_in.close();

    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;
    unsigned int dom_count = 0;
    double avg_relative_error = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    double relative_error = effectiveness::COD_prop_personalized_cross_precision(qp, &g, std::stod(topr), beta, visited, back_visited, method, hg_doms, hg_isdoms, hg_doms_greater, hg_isdoms_greater, dom_count, avg_time);
                    if(relative_error >= -0.1){
                        rule_count++;
                        avg_relative_error += relative_error;
                    }

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    double relative_error = effectiveness::COD_prop_personalized_cross_precision(qp, &g, std::stod(topr), beta, visited, back_visited, method, hg_doms, hg_isdoms, hg_doms_greater, hg_isdoms_greater, dom_count, avg_time);
                    if(relative_error > -0.1){
                        rule_count++;
                        avg_relative_error += relative_error;
                    }
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_relative_error /= rule_count;
    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~goodness:"<<avg_relative_error<<std::endl;
    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Scalability_hg_greater_f1(const std::string &choice, double topr, const std::string &method, const std::string &len){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat"+len);
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    unsigned int rule_count = 0;
    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool peer_found = effectiveness::COD_hg_global_greater_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                    if(peer_found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool peer_found = effectiveness::COD_hg_global_greater_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                    if(peer_found) rule_count++;
                    if(rule_count >= PATHCOUNT) break;

                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_hg_global_greater_f1(const std::string &choice, double topr, const std::string &method){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    unsigned int rule_count = 0;
    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
                int sub = stoi(rules_line.substr(i, temp_pos - i));

                if (sub == -1) state = 0; // next int indicates a variable rule
                else if (sub == -2) qp->EDirect.push_back(1); // ->
                else if (sub == -3) qp->EDirect.push_back(-1); // <-
                else if (sub == -4) { // pop
                    qp->EDirect.pop_back();
                    qp->ETypes.pop_back();
                    qp->NTypes.pop_back();
                } else if (sub == -5) state = 1; // next int indicates a instance rule
                else {
                    if (state == 1) { // instance rule
                        qp->instance = sub;
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        bool peer_found = effectiveness::COD_hg_global_greater_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                        if(peer_found) rule_count++;

                        qp->instance = -1;
                        state = -1;
                    } else if (state == 0) { // variable rule
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        bool peer_found = effectiveness::COD_hg_global_greater_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                        if(peer_found) rule_count++;

                        state = -1;
                    } else {
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                    }
                }
                i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_hg_global_f1_by_union(const std::string &choice, double topr, const std::string &method){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    unsigned int rule_count = 0;
    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool peer_found = effectiveness::COD_hg_global_f1_by_union(qp, &g, topr, visited, back_visited, method, avg_time);
                    if(peer_found) rule_count++;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    bool peer_found = effectiveness::COD_hg_global_f1_by_union(qp, &g, topr, visited, back_visited, method, avg_time);
                    if(peer_found) rule_count++;

                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }

    qp->clear();
    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_hg_global_f1(const std::string &choice, double topr, const std::string &method){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    unsigned int rule_count = 0;
    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
                int sub = stoi(rules_line.substr(i, temp_pos - i));

                if (sub == -1) state = 0; // next int indicates a variable rule
                else if (sub == -2) qp->EDirect.push_back(1); // ->
                else if (sub == -3) qp->EDirect.push_back(-1); // <-
                else if (sub == -4) { // pop
                    qp->EDirect.pop_back();
                    qp->ETypes.pop_back();
                    qp->NTypes.pop_back();
                } else if (sub == -5) state = 1; // next int indicates a instance rule
                else {
                    if (state == 1) { // instance rule
                        qp->instance = sub;
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        bool peer_found = effectiveness::COD_hg_global_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                        if(peer_found) rule_count++;

                        qp->instance = -1;
                        state = -1;
                    } else if (state == 0) { // variable rule
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        bool peer_found = effectiveness::COD_hg_global_f1(qp, &g, topr, visited, back_visited, method, avg_time);
                        if(peer_found) rule_count++;

                        state = -1;
                    } else {
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                    }
                }
                i = temp_pos;
        }
    }

    qp->clear();
    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_time /= rule_count;
    avg_time /= 1000000000;

    std::cout<<"~time_per_rule:"<<avg_time<<" s"<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}

void Effective_hg_stats(const std::string &choice){
    HeterGraph g(choice);
    std::string qnodes_path = (choice + "/qnodes_"+choice+".dat");
    std::ifstream qnodes_in;
    qnodes_in.open(qnodes_path);
    std::vector<unsigned int> qnodes;
    std::string qnode_line;
    while(getline(qnodes_in, qnode_line)){
        auto qnode = (unsigned int)(stoi(qnode_line));
        qnodes.push_back(qnode);
    }
    qnodes_in.close();

    std::string rules_path = (choice + "/cod-rules_"+choice+".limit");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    int query_count = -1;
    std::vector<double> res;
    for(unsigned int i=0;i < 4;i++) res.push_back(0.0);

    while(getline(rules_in, rules_line)){
        query_count++;
        unsigned int qnode = qnodes[query_count];
        int state = -1;
        // state == 0: next int represents a variable rule;
        // state == 1: next int represents a instance rule;

        std::string::size_type temp_pos;
        int size = rules_line.size();

        std::cout<<"@ qn:"<<qnode<<std::endl;

        std::vector<double> qn_res;
        for(unsigned int i=0;i<res.size();i++) qn_res.push_back(0.0);
        int qn_rule_count = 0;

        for(unsigned int i = 0; i<size;i++) {
            temp_pos = rules_line.find(' ', i);
            if (temp_pos < size) {
                int sub = stoi(rules_line.substr(i, temp_pos - i));

                if (sub == -1) state = 0; // next int indicates a variable rule
                else if (sub == -2) qp->EDirect.push_back(1); // ->
                else if (sub == -3) qp->EDirect.push_back(-1); // <-
                else if (sub == -4) { // pop
                    qp->EDirect.pop_back();
                    qp->ETypes.pop_back();
                    qp->NTypes.pop_back();
                } else if (sub == -5) state = 1; // next int indicates a instance rule
                else {
                    if (state == 1) { // instance rule
                        qp->instance = sub;
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        auto temp_res = effectiveness::COD_hg_statistics(qp, &g, qnode, visited, back_visited);
                        if(temp_res != nullptr) {
                            for (unsigned int z = 0; z < res.size(); z++) qn_res[z] += temp_res->at(z);
                            qn_rule_count++;
                            delete temp_res;
                        }
                        qp->instance = -1;
                        state = -1;
                    } else if (state == 0) { // variable rule
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                        while(visited->size()< qp->ETypes.size()+1){
                            visited->push_back(new std::vector<bool>(g.NT.size(), false));
                            back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        }
                        auto temp_res = effectiveness::COD_hg_statistics(qp, &g, qnode, visited, back_visited);
                        if(temp_res != nullptr) {
                            for (unsigned int z = 0; z < res.size(); z++) qn_res[z] += temp_res->at(z);
                            qn_rule_count++;
                            delete temp_res;
                        }
                        state = -1;
                    } else {
                        qp->ETypes.push_back(sub);
                        qp->NTypes.push_back(-1);
                    }
                }
                i = temp_pos;
            }
        }
        qp->clear();
        for(unsigned int i=0;i<res.size();i++) {
            qn_res[i] /= qn_rule_count;
            res[i] += qn_res[i];
        }
        std::cout<<"q_dens:"<<qn_res[0]<<std::endl;
        std::cout<<"q_d_same:"<<qn_res[1]<<std::endl;
        std::cout<<"q_h_same:"<<qn_res[2]<<std::endl;
        std::cout<<"q_|peer|:"<<qn_res[3]<<std::endl;
    }
    rules_in.close();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    for(unsigned int i=0;i<res.size();i++) res[i] /= qnodes.size();
    std::cout<<"~dens:"<<res[0]<<std::endl;
    std::cout<<"~d_same:"<<res[1]<<std::endl;
    std::cout<<"~h_same:"<<res[2]<<std::endl;
    std::cout<<"~|peer|:"<<res[3]<<std::endl;
}