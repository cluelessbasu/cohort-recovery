#ifndef ESCAPE_DECOMP_H_
#define ESCAPE_DECOMP_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/ClusterStructures.h"

#include <omp.h>
#include <stack>
#include <queue>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <cstring>
#include <unordered_set>
#include <iomanip>

using namespace Escape;
using namespace std;

const int NUM_THREADS = 12;

struct DiagnosticStats {
    long total_attempts = 0;       
    long any_nonzero_event = 0;    
    
    long count_only_m1 = 0;        
    long count_only_m2 = 0;        
    long count_both = 0;           
    
    long long sum_size_m1 = 0;
    long long sum_size_m2 = 0;
    
    long count_m1_nonzero = 0;
    long count_m2_nonzero = 0;

    long long sum_intersection = 0;
    long intersection_samples = 0;

    double sum_new_from_core = 0;
    long count_core_adds_any = 0;
};

void populateStats(CGraph* g, Cluster* clus)
{
    clus->nVertices = (clus->vertices).size(); 
    clus->nEdges = 0; 
    clus->cut = 0; 
    for (auto iter = begin(clus->vertices); iter != end(clus->vertices); iter++) 
    {
        VertexIdx u = *iter; 
        for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; i++) 
        {
            VertexIdx v = g->nbors[i]; 
            if ((clus->vertices).find(v) != (clus->vertices).end()) 
                clus->nEdges++; 
            else 
                clus->cut++;
        }
    }
    clus->nEdges = clus->nEdges/2; 
}

WeightedTriangleInfo commonNbr(CGraph *g)
{
   printf("Computing all triangle weighted info\n");

   WeightedTriangleInfo ret;   
   ret.total = 0;      
   ret.total_wgt = 0;  
   ret.perVertex = new Weight[g->nVertices+1];
   ret.perEdge = new Weight[g->nEdges+1]; 

   #pragma omp parallel for
   for (VertexIdx i=0; i < g->nVertices; ++i) ret.perVertex[i] = 0;
   #pragma omp parallel for
   for (EdgeIdx j=0; j < g->nEdges; ++j) ret.perEdge[j] = 0; 

   g->getPartnerMap(); 

   #pragma omp parallel
   {
        WeightedTriangleInfo local;
        local.total = 0;
        local.total_wgt = 0;
        local.perVertex = new Weight[g->nVertices+1]();
        local.perEdge = new Weight[g->nEdges+1]();

        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx u=0; u < g->nVertices; ++u) 
            for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; ++i)   
            {
                VertexIdx v = g->nbors[i]; 
                Count degu = g->offsets[u+1] - g->offsets[u]; 
                Count degv = g->offsets[v+1] - g->offsets[v]; 

                if(u > v) continue;

                VertexIdx lower = (degu < degv) ? u : v; 
                VertexIdx higher = (degu < degv)? v : u; 

                for (EdgeIdx j = g->offsets[lower]; j < g->offsets[lower+1]; ++j)
                {
                    VertexIdx w = g->nbors[j]; 
                    if (w <= u || w <= v) continue;
                    
                    EdgeIdx loc = g->getEdgeBinary(higher,w); 
                    if (loc != -1) 
                    {
                        Count degw = g->offsets[w+1] - g->offsets[w]; 
                        Weight wgt = 1.0/double(degu*degv*degw); 
                        
                        local.total++; 
                        local.total_wgt += wgt; 

                        local.perEdge[i] += wgt; 
                        local.perEdge[g->partnerMap[i]] += wgt; 
                        local.perEdge[j] += wgt; 
                        local.perEdge[g->partnerMap[j]] += wgt; 
                        local.perEdge[loc] += wgt; 
                        local.perEdge[g->partnerMap[loc]] += wgt; 
                        
                        local.perVertex[u] += wgt; 
                        local.perVertex[v] += wgt; 
                        local.perVertex[w] += wgt; 
                    }
                }
            }
        
        #pragma omp critical
        {
            ret.total += local.total;
            ret.total_wgt += local.total_wgt;
            for (VertexIdx x = 0; x < g->nVertices; x++)
                ret.perVertex[x] += local.perVertex[x];
            for (EdgeIdx e = 0; e < g->nEdges; e++)
                ret.perEdge[e] += local.perEdge[e];
        }

        delete[] local.perVertex;
        delete[] local.perEdge;
   }
   return ret;
}

char* multiCleaner(CGraph* g, WeightedTriangleInfo* triInfo, double eps)
{
    char* edgeStatus = new char[g->nEdges];
    std::fill(edgeStatus, edgeStatus + g->nEdges, 'Y');

    stack<Pair> toDelete;

    for (VertexIdx u = 0; u < g->nVertices; ++u) {
        Count degu = g->offsets[u+1] - g->offsets[u];
        for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; ++i) {
            VertexIdx v = g->nbors[i];
            if (u >= v) continue;

            Count degv = g->offsets[v+1] - g->offsets[v];
            double thresh = eps * (1.0 / (double)(degu * degv));

            if (triInfo->perEdge[i] < thresh) {
                toDelete.push({u,v});
                edgeStatus[i] = 'S';
                edgeStatus[g->partnerMap[i]] = 'S';
            }
        }
    }

    while (!toDelete.empty()) {
        Pair e = toDelete.top(); toDelete.pop();
        VertexIdx u = e.first, v = e.second;

        EdgeIdx uv = g->getEdgeBinary(u,v);
        if (uv == -1) continue;
        if (edgeStatus[uv] == 'N') continue;

        EdgeIdx vu = g->partnerMap[uv];
        edgeStatus[uv] = edgeStatus[vu] = 'N';

        Count degu = g->offsets[u+1] - g->offsets[u];
        Count degv = g->offsets[v+1] - g->offsets[v];

        VertexIdx low = (degu < degv ? u : v);
        VertexIdx high = (degu < degv ? v : u);

        for (EdgeIdx j = g->offsets[low]; j < g->offsets[low+1]; ++j) {
            if (edgeStatus[j] != 'Y') continue;
            VertexIdx w = g->nbors[j];
            if (w == u || w == v) continue;

            EdgeIdx hw = g->getEdgeBinary(high, w);
            if (hw == -1 || edgeStatus[hw] != 'Y') continue;

            Count degw = g->offsets[w+1] - g->offsets[w];
            double wgt = 1.0 / (double)(degu * degv * degw);

            triInfo->perEdge[j] -= wgt;
            triInfo->perEdge[g->partnerMap[j]] -= wgt;
            triInfo->perEdge[hw] -= wgt;
            triInfo->perEdge[g->partnerMap[hw]] -= wgt;

            double thresh_l = eps * (1.0 / (double)(min(degu,degv) * degw));
            if (triInfo->perEdge[j] < thresh_l && edgeStatus[j] == 'Y') {
                edgeStatus[j] = edgeStatus[g->partnerMap[j]] = 'S';
                toDelete.push({low,w});
            }

            double thresh_h = eps * (1.0 / (double)(max(degu,degv) * degw));
            if (triInfo->perEdge[hw] < thresh_h && edgeStatus[hw] == 'Y') {
                edgeStatus[hw] = edgeStatus[g->partnerMap[hw]] = 'S';
                toDelete.push({high,w});
            }
        }
    }
    return edgeStatus;
}

int runCoreOptimized(CGraph* g, 
                      const char* edgeStatus,
                      const std::vector<VertexIdx>& candidate_vertices, 
                      Cluster& target_cluster, 
                      std::vector<int>& degrees,           
                      std::vector<std::vector<int>>& adj,  
                      std::vector<bool>& removed,          
                      std::vector<int>& queue,             
                      std::vector<int>& global_to_local,
                      std::atomic<bool>* core_tracker,
                      int k_core_val = 2) 
{
    std::vector<VertexIdx> local_to_global;
    local_to_global.reserve(candidate_vertices.size());

    for (VertexIdx u : candidate_vertices) {
        global_to_local[u] = local_to_global.size();
        local_to_global.push_back(u);
    }

    int n = local_to_global.size();
    if (n == 0) {
        for(auto u : local_to_global) global_to_local[u] = -1; 
        return 0;
    }

    if (degrees.size() < n) degrees.resize(n);
    if (adj.size() < n) adj.resize(n);
    if (removed.size() < n) removed.resize(n);

    for(int i=0; i<n; ++i) {
        degrees[i] = 0;
        adj[i].clear();
        removed[i] = false;
    }

    for (int i = 0; i < n; i++) {
        VertexIdx u = local_to_global[i];
        for (EdgeIdx j = g->offsets[u]; j < g->offsets[u+1]; j++) {
            if (edgeStatus[j] != 'Y') continue;
            VertexIdx v = g->nbors[j];
            int neighbor_idx = global_to_local[v];
            
            if (neighbor_idx != -1 && neighbor_idx < n && local_to_global[neighbor_idx] == v) {
                adj[i].push_back(neighbor_idx);
                degrees[i]++;
            }
        }
    }

    queue.clear();
    int head = 0;

    for(int i=0; i<n; i++) {
        if(degrees[i] < k_core_val) {
            removed[i] = true;
            queue.push_back(i);
        }
    }

    while(head < queue.size()){
        int u_local = queue[head++];
        for(int v_local : adj[u_local]){
            if(!removed[v_local]){
                degrees[v_local]--;
                if(degrees[v_local] < k_core_val){
                    removed[v_local] = true;
                    queue.push_back(v_local);
                }
            }
        }
    }

    int added_count = 0;
    for(int i=0; i<n; i++){
        if(!removed[i]){
            if(target_cluster.vertices.insert(local_to_global[i]).second) {
                added_count++;
                if (core_tracker) {
                    core_tracker[local_to_global[i]].store(true, std::memory_order_relaxed);
                }
            }
        }
    }

    for(auto u : local_to_global) global_to_local[u] = -1;
    return added_count;
}

void checkCleaningIntegrity(CGraph* g, char* edgeStatus, WeightedTriangleInfo* triInfo, double eps)
{
    printf("\n--- VERIFYING CLEANING STATUS ---\n");
    long active_edges = 0;
    long asymmetric_count = 0;
    long violation_count = 0; 

    #pragma omp parallel for reduction(+:active_edges, asymmetric_count, violation_count)
    for (VertexIdx u = 0; u < g->nVertices; ++u) 
    {
        for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; ++i) 
        {
            VertexIdx v = g->nbors[i];
            EdgeIdx partner = g->partnerMap[i];
            if (edgeStatus[i] != edgeStatus[partner]) {
                asymmetric_count++;
            }
            if (edgeStatus[i] == 'Y') {
                active_edges++;
                Count degu = g->offsets[u+1] - g->offsets[u];
                Count degv = g->offsets[v+1] - g->offsets[v];
                double edgeWt = 1.0 / (double)(degu * degv);
                
                if (triInfo->perEdge[i] < (eps * 0.2 * edgeWt) - 1e-9) {
                    violation_count++;
                }
            }
        }
    }

    printf("Active Edges ('Y'): %ld\n", active_edges);
    printf("1. Asymmetry Errors: %ld\n", asymmetric_count / 2);
    if (asymmetric_count > 0) printf("   [FAIL] Inconsistent graph topology.\n");
    else printf("   [PASS] Graph topology is symmetric.\n");

    printf("2. Stability Violations: %ld\n", violation_count);
    if (violation_count > 0) printf("   [FAIL] Cleaning did not converge.\n");
    else printf("   [PASS] All edges satisfy condition.\n");
    printf("---------------------------------\n\n");
}

vector<Cluster> IndExtract(CGraph* g, double eps, std::map<VertexIdx, VertexIdx>& cluster_map)
{ 
    auto start = std::chrono::high_resolution_clock::now();
    vector<Cluster> decomposition; 

    for (VertexIdx i=0; i < g->nVertices; i++) {
        cluster_map[i] = -1; 
    }
    
    WeightedTriangleInfo triInfo = commonNbr(g); 
    auto stop = std::chrono::high_resolution_clock::now();
    cout << "Time for weights: " << chrono::duration<double>(stop - start).count() << "s \n";
    start = std::chrono::high_resolution_clock::now();

    char* edgeStatus = multiCleaner(g, &triInfo, eps);
    checkCleaningIntegrity(g, edgeStatus, &triInfo, eps); 

    stop = std::chrono::high_resolution_clock::now();
    cout << "Time for initial clean: " << chrono::duration<double>(stop - start).count() << "s \n";
    start = std::chrono::high_resolution_clock::now();

    std::vector<Pair> deg_info(g->nVertices);
    #pragma omp parallel for
    for (VertexIdx i=0; i < g->nVertices; i++) {
        deg_info[i].first = i;
        deg_info[i].second = g->offsets[i+1] - g->offsets[i];
    }
    std::sort(deg_info.begin(), deg_info.end(), [](const Pair& a, const Pair& b) {
        if (a.second != b.second) return a.second < b.second; 
        return a.first < b.first; 
    }); 
    
    std::vector<VertexIdx> mis_seeds;
    std::vector<bool> covered(g->nVertices, false);
    
    for(size_t i=0; i < g->nVertices; ++i) {
        VertexIdx u = deg_info[i].first;
        if (covered[u]) continue; 
        
        bool u_active = false;
        for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) {
            if (edgeStatus[j] == 'Y') { u_active = true; break; }
        }
        if (!u_active) continue;

        mis_seeds.push_back(u);
        covered[u] = true;
        
        Count degu = g->offsets[u+1] - g->offsets[u];
        for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) {
            if (edgeStatus[j] == 'Y') {
                VertexIdx v = g->nbors[j];
                Count degv = g->offsets[v+1] - g->offsets[v];
                if ((double)degv <= (double)degu/eps) {
                    covered[v] = true;
                }
            }
        }
    }

    stop = std::chrono::high_resolution_clock::now();
    cout << "MIS Construction Complete. Seeds found: " << mis_seeds.size() << "\n";
    cout << "Time for MIS: " << chrono::duration<double>(stop - start).count() << "s \n";
    start = std::chrono::high_resolution_clock::now();

    std::atomic<bool>* core_tracker = new std::atomic<bool>[g->nVertices];
    for(size_t i=0; i<g->nVertices; ++i) std::atomic_init(&core_tracker[i], false);

    DiagnosticStats global_stats;

    #pragma omp parallel 
    {
        std::vector<double> triwgt(g->nVertices, 0.0);
        std::vector<bool> nbdBitMap(g->nVertices, false);
        std::vector<bool> isTwoHop(g->nVertices, false);
        
        std::vector<VertexIdx> two_hop_vec; two_hop_vec.reserve(1024);
        std::vector<wgtPair> candidates; candidates.reserve(1024);
        std::vector<std::pair<double, std::pair<VertexIdx, VertexIdx>>> edgeCandidates; edgeCandidates.reserve(1024);
        std::vector<VertexIdx> core_candidates; core_candidates.reserve(1024);
        std::vector<VertexIdx> m1_added, m2_added;

        std::vector<int> rc_degrees;
        std::vector<std::vector<int>> rc_adj;
        std::vector<bool> rc_removed;
        std::vector<int> rc_queue;
        std::vector<int> rc_map(g->nVertices, -1);

        std::vector<Cluster> local_decomposition; 
        DiagnosticStats local_stats;

        #pragma omp for schedule(dynamic, 64)
        for (size_t i=0; i < mis_seeds.size(); i++) 
        {
            VertexIdx u = mis_seeds[i];
            
            two_hop_vec.clear();
            candidates.clear();
            edgeCandidates.clear();
            m1_added.clear();
            m2_added.clear();

            Count degu = g->offsets[u+1] - g->offsets[u];
            Cluster next_clus; 
            next_clus.vertices.insert(u); 
            double internalWgt = 0; 

            for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) {
                VertexIdx v = g->nbors[j];
                if (edgeStatus[j] != 'Y') continue;
                next_clus.vertices.insert(v);
                nbdBitMap[v] = true; 
            }

            for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) {
                VertexIdx v = g->nbors[j];
                if (edgeStatus[j] != 'Y') continue;
                Count degv = g->offsets[v+1] - g->offsets[v];
                if ((double)degv > (double)degu/eps) continue;
                for (EdgeIdx k=j; k < g->offsets[u+1]; ++k) {
                    VertexIdx w = g->nbors[k];
                    if (edgeStatus[k] != 'Y') continue; 
                    Count degw = g->offsets[w+1] - g->offsets[w];
                    if ((double)degw > (double)degu/eps) continue;
                    EdgeIdx loc = g->getEdgeBinary(u,w); 
                    if (loc == -1 || edgeStatus[loc] != 'Y') continue; 
                    internalWgt += 1.0/(double)(degu*degv*degw);
                    VertexIdx lower = (degv < degw) ? v : w;
                    VertexIdx higher = (degv < degw) ? w : v;
                    for (EdgeIdx ell = g->offsets[lower]; ell < g->offsets[lower+1]; ell++) {
                        if(edgeStatus[ell] != 'Y') continue; 
                        VertexIdx x = g->nbors[ell];
                        if(x == lower || x == higher || x == u) continue;
                        EdgeIdx loc_hix = g->getEdgeBinary(higher,x); 
                        if(loc_hix == -1 || edgeStatus[loc_hix] != 'Y') continue; 
                        Count degx = g->offsets[x+1] - g->offsets[x];
                        double wgt = 1.0/(double)(degx*degv*degw);
                        if (nbdBitMap[x]) {
                            if (x > v && x > w) internalWgt += wgt; 
                            continue; 
                        }
                        if (!isTwoHop[x]) {
                            two_hop_vec.push_back(x);
                            isTwoHop[x] = true;
                        }
                        triwgt[x] += wgt; 
                    }
                }
            }

            double potentialWgt = 0;
            for (VertexIdx curr : two_hop_vec) {
                candidates.push_back({curr, triwgt[curr]});
                potentialWgt += triwgt[curr];
                triwgt[curr] = 0; 
            }
            std::sort(candidates.begin(), candidates.end(), [](const wgtPair& a, const wgtPair& b) {
                return a.wgt > b.wgt;
            });

            double ratio = internalWgt/(internalWgt + potentialWgt);
            int maxind = -1;
            double num_sum = 0, den_sum = 0;
            for (int idx = 0; idx < candidates.size(); idx++) {
                num_sum += candidates[idx].wgt;
                den_sum += triInfo.perVertex[candidates[idx].vertex];
                double new_ratio = (internalWgt + num_sum)/(internalWgt + potentialWgt + den_sum);
                if (new_ratio > ratio) {
                    ratio = new_ratio;
                    maxind = idx;
                }
            }
            for (int idx=0; idx <= maxind; idx++) {
                VertexIdx next_vert = candidates[idx].vertex;
                next_clus.vertices.insert(next_vert);
                m1_added.push_back(next_vert);
            }
            std::sort(m1_added.begin(), m1_added.end()); 
            
            for (auto v : two_hop_vec) {
                for (EdgeIdx j = g->offsets[v]; j < g->offsets[v+1]; ++j) {
                    VertexIdx w = g->nbors[j];
                    if (edgeStatus[j] != 'Y') continue;
                    if (w <= v) continue; 
                    if (!isTwoHop[w]) continue; 
                    double edgeScore = 0.0;
                    EdgeIdx idx_v = g->offsets[v];
                    EdgeIdx idx_w = g->offsets[w];
                    EdgeIdx end_v = g->offsets[v+1];
                    EdgeIdx end_w = g->offsets[w+1];
                    while(idx_v < end_v && idx_w < end_w) {
                        VertexIdx nv = g->nbors[idx_v];
                        VertexIdx nw = g->nbors[idx_w];
                        if (nv == nw) {
                            if (nbdBitMap[nv]) {
                                Count deg_z = g->offsets[nv+1] - g->offsets[nv]; 
                                Count deg_v_curr = g->offsets[v+1] - g->offsets[v];
                                Count deg_w_curr = g->offsets[w+1] - g->offsets[w];
                                edgeScore += 1.0 / (double)(deg_v_curr * deg_w_curr * deg_z);
                            }
                            idx_v++; idx_w++;
                        } else if (nv < nw) idx_v++;
                        else idx_w++;
                    }
                    if (edgeScore > 0) edgeCandidates.push_back({edgeScore, {v, w}});
                }
            }
            std::sort(edgeCandidates.begin(), edgeCandidates.end(), 
                      std::greater<std::pair<double, std::pair<VertexIdx, VertexIdx>>>());

            int maxEdgeInd = -1;
            double edge_num_sum = 0;
            double edge_den_sum = 0;
            std::unordered_set<VertexIdx> distinct_in_prefix; 
            for (int i = 0; i < edgeCandidates.size(); i++) {
                double score = edgeCandidates[i].first;
                VertexIdx v = edgeCandidates[i].second.first;
                VertexIdx w = edgeCandidates[i].second.second;
                edge_num_sum += score;
                double cost_v = 0, cost_w = 0;
                bool v_in_cluster = (next_clus.vertices.find(v) != next_clus.vertices.end());
                bool w_in_cluster = (next_clus.vertices.find(w) != next_clus.vertices.end());
                bool v_in_prefix = (distinct_in_prefix.find(v) != distinct_in_prefix.end());
                bool w_in_prefix = (distinct_in_prefix.find(w) != distinct_in_prefix.end());
                if (!v_in_cluster && !v_in_prefix) cost_v = triInfo.perVertex[v];
                if (!w_in_cluster && !w_in_prefix) cost_w = triInfo.perVertex[w];
                edge_den_sum += (cost_v + cost_w);
                double new_ratio = (internalWgt + num_sum + edge_num_sum) / (internalWgt + potentialWgt + den_sum + edge_den_sum);
                if (new_ratio > ratio) {
                    ratio = new_ratio;
                    maxEdgeInd = i;
                }
                if (!v_in_cluster) distinct_in_prefix.insert(v);
                if (!w_in_cluster) distinct_in_prefix.insert(w);
            }

            if (maxEdgeInd != -1) {
                core_candidates.clear();
                for(int i=0; i <= maxEdgeInd; ++i) {
                     core_candidates.push_back(edgeCandidates[i].second.first);
                     core_candidates.push_back(edgeCandidates[i].second.second);
                }
                std::sort(core_candidates.begin(), core_candidates.end());
                core_candidates.erase(std::unique(core_candidates.begin(), core_candidates.end()), core_candidates.end());
                m2_added = core_candidates;
                int64_t pre_core_size = next_clus.vertices.size();
                
                int new_from_core = runCoreOptimized(g, edgeStatus, core_candidates, next_clus, rc_degrees, rc_adj, rc_removed, rc_queue, rc_map, core_tracker);
                
                if (new_from_core > 0) {
                    local_stats.sum_new_from_core += double(new_from_core)/double(pre_core_size);
                    local_stats.count_core_adds_any++;
                }
            }

            local_stats.total_attempts++;
            bool m1_nonzero = !m1_added.empty();
            bool m2_nonzero = !m2_added.empty();
            if (m1_nonzero) {
                local_stats.sum_size_m1 += m1_added.size();
                local_stats.count_m1_nonzero++;
            }
            if (m2_nonzero) {
                local_stats.sum_size_m2 += m2_added.size();
                local_stats.count_m2_nonzero++;
            }
            if (m1_nonzero || m2_nonzero) {
                local_stats.any_nonzero_event++;
                if (m1_nonzero && !m2_nonzero) local_stats.count_only_m1++;
                else if (!m1_nonzero && m2_nonzero) local_stats.count_only_m2++;
                else if (m1_nonzero && m2_nonzero) {
                    local_stats.count_both++;
                    std::vector<VertexIdx> intersection;
                    std::set_intersection(m1_added.begin(), m1_added.end(),
                                          m2_added.begin(), m2_added.end(),
                                          std::back_inserter(intersection));
                    local_stats.sum_intersection += intersection.size();
                    local_stats.intersection_samples++;
                }
            }

            for(auto v : two_hop_vec) isTwoHop[v] = false;
            for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) nbdBitMap[g->nbors[j]] = false;

            if (next_clus.vertices.size() >= 2) 
            {
                 populateStats(g, &next_clus); 
                 local_decomposition.push_back(next_clus);
            }
        } 

        #pragma omp critical
        {
            global_stats.total_attempts += local_stats.total_attempts;
            global_stats.any_nonzero_event += local_stats.any_nonzero_event;
            global_stats.count_only_m1 += local_stats.count_only_m1;
            global_stats.count_only_m2 += local_stats.count_only_m2;
            global_stats.count_both += local_stats.count_both;
            global_stats.sum_size_m1 += local_stats.sum_size_m1;
            global_stats.sum_size_m2 += local_stats.sum_size_m2;
            global_stats.count_m1_nonzero += local_stats.count_m1_nonzero;
            global_stats.count_m2_nonzero += local_stats.count_m2_nonzero;
            global_stats.sum_intersection += local_stats.sum_intersection;
            global_stats.intersection_samples += local_stats.intersection_samples;
            global_stats.sum_new_from_core += local_stats.sum_new_from_core;
            global_stats.count_core_adds_any += local_stats.count_core_adds_any;

            int start_idx = decomposition.size();
            decomposition.insert(decomposition.end(), local_decomposition.begin(), local_decomposition.end());
            for(size_t c=0; c<local_decomposition.size(); ++c) {
                for(auto v : local_decomposition[c].vertices) {
                    cluster_map[v] = start_idx + c;
                }
            }
        }
    } 

    stop = std::chrono::high_resolution_clock::now();
    cout << "Time for clusters: " << chrono::duration<double>(stop - start).count() << "s \n";

    cout << "\n--- MIS-Based Extraction Diagnostics (Non-Zero Cases) ---\n";
    cout << "Total Seeds (MIS Size): " << global_stats.total_attempts << "\n";
    cout << "Seeds with ANY 2-hop contribution: " << global_stats.any_nonzero_event << "\n\n";

    long any = global_stats.any_nonzero_event > 0 ? global_stats.any_nonzero_event : 1;
    cout << "1. Contribution Type (Percentage of Non-Zero events):\n";
    cout << "   Vertex-Based ONLY: " << global_stats.count_only_m1 << " (" << (100.0 * global_stats.count_only_m1 / any) << "%)\n";
    cout << "   Edge-Based ONLY:   " << global_stats.count_only_m2 << " (" << (100.0 * global_stats.count_only_m2 / any) << "%)\n";
    cout << "   Both Methods:      " << global_stats.count_both << " (" << (100.0 * global_stats.count_both / any) << "%)\n\n";

    double avg_m1 = (double)global_stats.sum_size_m1 / (global_stats.count_m1_nonzero > 0 ? global_stats.count_m1_nonzero : 1);
    double avg_m2 = (double)global_stats.sum_size_m2 / (global_stats.count_m2_nonzero > 0 ? global_stats.count_m2_nonzero : 1);
    cout << "2. Average Size (Only counting non-zero cases):\n";
    cout << "   Method 1 (Vertex): " << avg_m1 << "\n";
    cout << "   Method 2 (Edge):   " << avg_m2 << "\n\n";
    
    double avg_int = (double)global_stats.sum_intersection / (global_stats.intersection_samples > 0 ? global_stats.intersection_samples : 1);
    cout << "3. Average Intersection (When both > 0): " << avg_int << "\n\n";
    
    double avg_core_add = (double)global_stats.sum_new_from_core / (global_stats.count_core_adds_any > 0 ? global_stats.count_core_adds_any : 1);
    cout << "4. Core Value-Add:\n";
    cout << "   Average NEW vertices added by Core (when it adds any): " << fixed << setprecision(2) << avg_core_add * 100 << "%\n";
    cout << "--------------------------------\n\n";

    long long total_cluster_size_sum = 0;
    std::vector<int> frequency(g->nVertices, 0);
    std::vector<bool> unique_clustered(g->nVertices, false);
    long unique_clustered_count = 0;

    for (const auto& clus : decomposition) {
        total_cluster_size_sum += clus.vertices.size();
        for (auto v : clus.vertices) {
            frequency[v]++;
            if (!unique_clustered[v]) {
                unique_clustered[v] = true;
                unique_clustered_count++;
            }
        }
    }

    long active_vertices = 0;
    #pragma omp parallel for reduction(+:active_vertices)
    for (VertexIdx u = 0; u < g->nVertices; ++u) {
        bool active = false;
        for (EdgeIdx j = g->offsets[u]; j < g->offsets[u+1]; ++j) {
            if (edgeStatus[j] == 'Y') { active = true; break; }
        }
        if (active) active_vertices++;
    }

    long core_overlap_count = 0;
    long total_core_added = 0;
    for (VertexIdx v=0; v < g->nVertices; ++v) {
        if (core_tracker[v].load()) {
            total_core_added++;
            if (frequency[v] > 1) core_overlap_count++;
        }
    }

    cout << "--- Decomposition Integrity Statistics ---\n";
    cout << "Sum of Cluster Sizes: " << total_cluster_size_sum << "\n";
    cout << "True Unique Clustered: " << unique_clustered_count << "\n";
    cout << "Active Vertices: " << active_vertices << "\n";
    cout << "Overlap Count: " << (total_cluster_size_sum - unique_clustered_count) << "\n";
    
    if (unique_clustered_count > active_vertices) {
        cout << "WARNING: GHOST NODES DETECTED! (Diff: " << (unique_clustered_count - active_vertices) << ")\n";
    } else {
        cout << "Sanity Check PASSED: Unique clustered vertices <= Active vertices.\n";
    }
    cout << "Unclustered Active Vertices: " << (active_vertices - unique_clustered_count) << "\n";
    
    cout << "\n--- CORE OVERLAP ANALYSIS ---\n";
    cout << "Unique vertices added by Core steps: " << total_core_added << "\n";
    cout << "Core vertices found in multiple clusters: " << core_overlap_count << "\n";
    if (total_core_added > 0)
        cout << "Percentage of Core vertices in overlaps: " << fixed << setprecision(2) << (100.0 * core_overlap_count / total_core_added) << "%\n";
    
    cout << "------------------------------------------\n";

    delete[] core_tracker;
    return decomposition;
}

vector<Cluster> connectedComp(CGraph* g)
{
    vector<Cluster> decomposition;
    bool* visited = new bool[g->nVertices]; 
    for (VertexIdx u=0; u < g->nVertices; u++) 
        visited[u] = false;

    queue<VertexIdx> q; 

    for (VertexIdx u=0; u < g->nVertices; u++) 
    {
        if (visited[u]) continue; 
        Cluster next_clus; 
        next_clus.vertices.insert(u); 
        visited[u] = true; 
        q.push(u); 

        while(!q.empty()) 
        {
            VertexIdx next = q.front(); 
            q.pop(); 
            for (EdgeIdx j = g->offsets[next]; j < g->offsets[next+1]; j++) 
            {
                VertexIdx nbr = g->nbors[j]; 
                if (!visited[nbr]) 
                {
                    q.push(nbr); 
                    visited[nbr] = true; 
                    next_clus.vertices.insert(nbr); 
                }
            }
        }
        populateStats(g,&next_clus); 
        decomposition.push_back(next_clus); 
    }
    delete[] visited;
    return decomposition;
}

#endif