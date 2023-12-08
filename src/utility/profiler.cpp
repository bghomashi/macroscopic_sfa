
#include "utility/profiler.h"
#include "utility/timer.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// this object to use as a log record.
struct ProfileObject {
    std::string name;           // name of this section of code
    double cumulative_time;     // time we have spent in this section of code
    int calls;                  // how many calls
    Timer timer;
};

static Profiler::Ptr_t s_profiler(new Profiler());
static std::unordered_map<std::string, ProfileObject> s_profiles = std::unordered_map<std::string, ProfileObject>();


void Profiler::Push(const std::string& name) {
    auto& profile = s_profiles[name];
    if (profile.name == "") {                   // set up first time
        profile.name = name;
        profile.cumulative_time = 0;
        profile.calls = 0;
    }
    profile.calls++;
    profile.timer.Reset();
}
void Profiler::Pop(const std::string& name) {
    auto& profile = s_profiles[name];
    profile.cumulative_time += profile.timer.Elapsed();
}
void Profiler::Reset() {
    s_profiles.clear();
}

bool Profiler::PrintTo(const std::string& filename) {
    std::vector<ProfileObject> profiles; profiles.reserve(s_profiles.size());

    for (const auto& p : s_profiles)
        profiles.push_back(p.second);
    
    std::sort(profiles.begin(), profiles.end(), [](const ProfileObject& a, const ProfileObject& b) {
        return a.cumulative_time < b.cumulative_time;
    });
    
    std::ofstream file(filename, std::ios_base::app);
    if (!file.is_open()) return false;

    file    << std::setw(32) << "Name" 
            << std::setw(10) << "Time [s]" 
            << std::setw(16) << "Calls" 
            << std::setw(16) << "[s]/Call" << std::endl;
    file << "-----------------------------------------------------------------" << std::endl;
    for (auto& p : profiles)
        file    << std::setw(32) << p.name 
                << std::setw(10) << p.cumulative_time 
                << std::setw(16) << p.calls
                << std::setw(16) << p.cumulative_time / p.calls << std::endl;
    file << "-----------------------------------------------------------------" << std::endl << std::endl;
    return true;
}
void Profiler::Print() {
    std::vector<ProfileObject> profiles; profiles.reserve(s_profiles.size());

    for (const auto& p : s_profiles)
        profiles.push_back(p.second);
    
    std::sort(profiles.begin(), profiles.end(), [](const ProfileObject& a, const ProfileObject& b) {
        return a.cumulative_time < b.cumulative_time;
    });
    std::cout   << std::setw(32) << "Name" 
                << std::setw(12) << "Time[s]" 
                << std::setw(16) << "Calls" 
                << std::setw(16) << "[s]/Call" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    for (auto& p : profiles)
        std::cout   << std::setw(32) << p.name 
                    << std::setw(12) << p.cumulative_time 
                    << std::setw(16) << p.calls
                    << std::setw(16) << p.cumulative_time / p.calls << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl << std::endl;
}


void Profile::Push(const std::string& name) {
    s_profiler->Push(name);
}
void Profile::Pop(const std::string& name) {
    s_profiler->Pop(name);
}

bool Profile::PrintTo(const std::string& filename) {
    return s_profiler->PrintTo(filename);
}
void Profile::Print() {
    s_profiler->Print();
}
void Profile::SetProfiler(Profiler* profiler) {
    s_profiler.reset(profiler);
}
void Profile::Reset() {
    s_profiler->Reset();
}

