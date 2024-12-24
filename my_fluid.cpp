
#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <ranges>
#include <cassert>
#include <random>
#include <algorithm>
#include <variant>
#include <tuple>
#include <charconv>
#include <fstream>
#include "FixedClass.h"
#include "Double.h"
#include "Float.h"
#ifdef FLOAT
#error "FLOAT is already defined"
#endif
#ifdef DOUBLE
#error "DOUBLE is already defined"
#endif
#ifdef FIXED
#error "FIXED is already defined"
#endif
#ifdef FAST_FIXED
#error "FAST_FIXED is already defined"
#endif
#ifdef S
#error "S is defined"
#endif

#define FLOAT float
#define DOUBLE double
#define FAST_FIXED(N, K) types::FastFixed<N, K>
#define FIXED(N, K)      types::Fixed<N, K>
#define S(N, M) types::sz{.rows = N, .columns = M, }
#define STRINGIFY(expr) #expr
inline constexpr std::string_view kFloatTypeName = STRINGIFY(FLOAT);
inline constexpr std::string_view kDoubleTypeName = STRINGIFY(DOUBLE);
inline constexpr std::string_view FFpref = "FAST_FIXED(";
inline constexpr std::string_view FFsuf = ")";
inline constexpr std::string_view Fpref = "FIXED(";
inline constexpr std::string_view Fsuf = ")";
constexpr size_t T = 1'000'000;
#undef STRINGIFY
namespace types {
    //Начинаю simulation

    struct sz {
        std::size_t rows{}, columns{};
    };

    struct Durex {
        std::size_t n, m;
        std::vector<std::vector<char>> field;
        double water_rho, air_rho, g;

        void resize() {
            field.resize(n);
            for(auto& str : field) {
                str.resize(m + 1);
            }
        }
    };

    template <class... ZipPack>
    struct List;

    template<std::size_t ind, class CurType, class... LeftTypes>
    struct Get_by_ind {
        using type = typename Get_by_ind<ind - 1, LeftTypes...>::type;
    };

    template<class CurType, class... LeftTypes>
    struct Get_by_ind<0, CurType, LeftTypes...> {
        using type = CurType;
    };


    template <std::size_t ind, class... LeftTypes>
    using Get_by_ind_t = typename Get_by_ind<ind, LeftTypes...>::type;

    constexpr std::size_t dyn_flag = std::numeric_limits<std::size_t>::max();

    template <class PPType, class VType, class VFlowType,
            std::size_t Rows = dyn_flag, std::size_t Columns = dyn_flag>
    class Simulator {
        static constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
        static constexpr bool static_flag = (Rows != dyn_flag) && (Columns != dyn_flag);
        //Для поля в котором на 1 больше столбик
        template <class EType>
         using StaticField = std::array<std::array<EType, Columns + 1>, Rows>;
        template <class EType>
         using DynamicField = std::vector<std::vector<EType>>;

        template <class EType>
        using FieldType = std::conditional_t<static_flag, StaticField<EType>, DynamicField<EType>>;

        //для обычного поля
        template <class EType>
        using StaticMatrix = std::array<std::array<EType, Columns>, Rows>;
        template <class EType>
        using DynamicMatrix = std::vector<std::vector<EType>>;

        template <class EType>
        using MatrixType = std::conditional_t<static_flag, StaticMatrix<EType>, DynamicMatrix<EType>>;

        template <class EType>
        struct VectorField {
            using CurMatrix = MatrixType<std::array<EType, deltas.size()>>;
            using CurType = CurMatrix::value_type;

            CurMatrix v{};
            VectorField() requires(static_flag) = default;

            VectorField(std::size_t cur_sz, const CurType& str) requires (!static_flag) :
            v{cur_sz, str} {}
            EType &add(int x, int y, int dx, int dy,const EType& dv) {
                return get(x, y, dx, dy) += dv;
            }
            EType &get(int x, int y, int dx, int dy) {
                size_t ind= std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
                assert(ind < deltas.size());
                return v[x][y][ind];
            }
        };
    public:
        using Field = FieldType<char>;
        using PPField = MatrixType<PPType>;
        using VField = VectorField<VType>;
        using VFlowField = VectorField<VFlowType>;
        using LastUseField = FieldType<int>;
        using dirsField = FieldType<int>;

        explicit constexpr Simulator(const Durex& cur) noexcept requires (static_flag) {
            for(std::size_t ind = 0; auto& str : field) {
                const std::vector<char>& cur_str = cur.field[ind];
                assert(cur_str.size() == str.size());
                std::ranges::copy(cur_str, str.begin());
            }
        }
        explicit constexpr Simulator(const Durex& cur)
        requires(!static_flag)
                : Simulator(cur, cur.n, cur.m - 1) {}
        [[nodiscard]] constexpr std::size_t get_n() const noexcept {
            return field.size();
        }
        [[nodiscard]] constexpr std::size_t get_m() const noexcept {
            return field.front().size() - 1;
        }
        std::tuple<PPType, bool, std::pair<int, int>> propagate_flow(int x, int y, PPType lim) {
            last_use_field[x][y] = UT - 1;
            auto ret = PPType(0);
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use_field[nx][ny] < UT) {
                    auto cap = v_field.get(x, y, dx, dy);
                    auto flow = v_flow_field.get(x, y, dx, dy);
                    if (VType(flow) == cap) {
                        continue;
                    }
                    // assert(v >= velocity_flow.get(x, y, dx, dy));
                    auto vp = types::min(VType(lim), cap - VType(flow));
                    if (last_use_field[nx][ny] == UT - 1) {
                        v_flow_field.add(x, y, dx, dy, VFlowType(vp));
                        last_use_field[x][y] = UT;
                        // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                        return {PPType(vp), 1, {nx, ny}};
                    }
                    auto [t, prop, end] = propagate_flow(nx, ny, PPType(vp));
                    ret += t;
                    if (prop) {
                        v_flow_field.add(x, y, dx, dy, VFlowType(t));
                        last_use_field[x][y] = UT;
                        // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                        return {t, prop && end != std::pair(x, y), end};
                    }
                }
            }
            last_use_field[x][y] = UT;
            return {ret, 0, {0, 0}};
        }

        void propagate_stop(int x, int y, bool force = false) {
            if (!force) {
                bool stop = true;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use_field[nx][ny] < UT - 1 && v_field.get(x, y, dx, dy) > VType(0)) {
                        stop = false;
                        break;
                    }
                }
                if (!stop) {
                    return;
                }
            }
            last_use_field[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use_field[nx][ny] == UT || v_field.get(x, y, dx, dy) > VType(0)) {
                    continue;
                }
                propagate_stop(nx, ny);
            }
        }

        PPType move_prob(int x, int y) {
            auto sum = PPType(0);
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use_field[nx][ny] == UT) {
                    continue;
                }
                auto v = v_field.get(x, y, dx, dy);
                if (v < VType(0)) {
                    continue;
                }
                sum += PPType(v);
            }
            return sum;
        }

        struct ParticleParams {
            char type;
            PPType cur_p;
            std::array<VType, deltas.size()> v;

            void swap_with(Simulator& simulator, int x, int y) {
                char field_cell = simulator.field[x][y];
                simulator.field[x][y] = type;
                type = field_cell;
                PPType p_copy = simulator.p[x][y];
                simulator.p[x][y] = cur_p;
                cur_p = p_copy;
                swap(simulator.v_field.v[x][y], v);
            }
        };

        double random01() {
            return (rnd() & ((1 << 16) - 1));
        }

        bool propagate_move(int x, int y, bool is_first) {
            last_use_field[x][y] = UT - is_first;
            bool ret = false;
            int nx = -1, ny = -1;
            do {
                std::array<PPType, deltas.size()> tres{};
                auto sum = PPType(0);
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use_field[nx][ny] == UT) {
                        tres[i] = sum;
                        continue;
                    }
                    auto v = v_field.get(x, y, dx, dy);
                    if (v < VType(0)) {
                        tres[i] = sum;
                        continue;
                    }
                    sum += PPType(v);
                    tres[i] = sum;
                }

                if (sum == PPType(0)) {
                    break;
                }

                PPType pp = PPType(random01()) * sum;
                size_t d = std::ranges::upper_bound(tres, pp) - tres.begin();

                auto [dx, dy] = deltas[d];
                nx = x + dx;
                ny = y + dy;
                assert(v_field.get(x, y, dx, dy) > VType(0) && field[nx][ny] != '#' && last_use_field[nx][ny] < UT);

                ret = (last_use_field[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
            } while (!ret);
            last_use_field[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use_field[nx][ny] < UT - 1 && v_field.get(x, y, dx, dy) < VType(0)) {
                    propagate_stop(nx, ny);
                }
            }
            if (ret) {
                if (!is_first) {
                    ParticleParams pp{};
                    pp.swap_with(*this, x, y);
                    pp.swap_with(*this, nx, ny);
                    pp.swap_with(*this, x, y);
                }
            }
            return ret;
        }

        void start(const Durex& cur) {
            rho[' '] = PPType(0.01);
            rho['.'] = PPType(1000);
            VType g(0.1);
            size_t N  = get_n();
            size_t M = get_m();
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                }
            }

            for (size_t i = 0; i < T; ++i) {
                PPType total_delta_p(0);
                // Apply external forces
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        if (field[x + 1][y] != '#')
                            v_field.add(x, y, 1, 0, g);
                    }
                }

                // Apply forces from p
                std::copy(p.begin(), p.end(), p_old.begin());
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            int nx = x + dx, ny = y + dy;
                            if (field[nx][ny] != '#' && p_old[nx][ny] < p_old[x][y]) {
                                auto delta_p = p_old[x][y] - p_old[nx][ny];
                                auto force = delta_p;
                                auto &contr = v_field.get(nx, ny, -dx, -dy);
                                if (contr * VType(rho[(int) field[nx][ny]]) >= VType(force)) {
                                    contr -= VType(force / rho[(int) field[nx][ny]]);
                                    continue;
                                }
                                force -= PPType(contr) * rho[(int) field[nx][ny]];
                                contr = VType(0);
                                v_field.add(x, y, dx, dy, (VType(force) / VType(rho[(int) field[x][y]])));
                                p[x][y] -= force / PPType(dirs[x][y]);
                                total_delta_p -= force / PPType(dirs[x][y]);
                            }
                        }
                    }
                }

                // Make flow from velocities
                bool prop = false;
                do {
                    UT += 2;
                    prop = false;
                    for (size_t x = 0; x < N; ++x) {
                        for (size_t y = 0; y < M; ++y) {
                            if (field[x][y] != '#' && last_use_field[x][y] != UT) {
                                auto [t, local_prop, _] = propagate_flow(x, y, PPType(1));
                                if (t > PPType(0)) {
                                    prop = 1;
                                }
                            }
                        }
                    }
                } while (prop);

                // Recalculate p with kinetic energy
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            auto old_v = v_field.get(x, y, dx, dy);
                            auto new_v = v_flow_field.get(x, y, dx, dy);
                            if (old_v > VType(0)) {
                                assert(VType(new_v) <= old_v);
                                v_field.get(x, y, dx, dy) = VType(new_v);
                                auto force = (old_v - VType(new_v)) * VType(rho[(int) field[x][y]]);
                                if (field[x][y] == '.')
                                    force *= VType(0.8);
                                if (field[x + dx][y + dy] == '#') {
                                    p[x][y] += PPType(force / VType(dirs[x][y]));
                                    total_delta_p += PPType(force / VType(dirs[x][y]));
                                } else {
                                    p[x + dx][y + dy] += PPType(force / VType(dirs[x + dx][y + dy]));
                                    total_delta_p += PPType(force / VType(dirs[x + dx][y + dy]));
                                }
                            }
                        }
                    }
                }

                UT += 2;
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use_field[x][y] != UT) {
                            if (PPType(random01()) < move_prob(x, y)) {
                                prop = true;
                                propagate_move(x, y, true);
                            } else {
                                propagate_stop(x, y, true);
                            }
                        }
                    }
                }
                    std::cout << "Tick " << i << ":\n";
                    for (size_t x = 0; x < N; ++x) {
                        for (size_t y = 0; y < M; ++y) {
                            std::cout << field[x][y];
                        }
                        std::cout << '\n';
                    }
            }
        }
    private:
        explicit constexpr Simulator(const Durex& cur, std::size_t rows, std::size_t columns) requires(!static_flag):
                field{cur.field},
                p{rows, typename PPField::value_type(columns)},
                p_old{rows, typename PPField::value_type(columns)},
                v_field{rows, typename VField::CurType(columns)},
                v_flow_field{rows, typename VFlowField::CurType(columns)},
                last_use_field{rows, typename LastUseField::value_type(columns)},
                rnd(1337),
                dirs{rows, typename dirsField::value_type(columns)}{}
        Field field{};
        PPField p, p_old{};
        VField v_field{};
        VFlowField v_flow_field{};
        LastUseField last_use_field{};
        dirsField dirs{};
        int UT = 0;
        std::array <PPType, 256> rho{};
        std::mt19937 rnd;
    };

    template <class PPType, class VType, class VFlowType,
            std::size_t Rows = dyn_flag, std::size_t Columns = dyn_flag>
    void start_on_field(const Durex& cur) {
        Simulator<PPType, VType, VFlowType, Rows, Columns> a{cur};
        a.start(cur);
    }

    template <class PPType, class VType, class VFlowType,
            sz CurSize, sz... LeftSizes>
    void select_size_and_start(const Durex& cur) {
        if(CurSize.rows == cur.field.size() && CurSize.columns == cur.field.front().size() - 1) {
            start_on_field<PPType, VType, VFlowType, CurSize.rows, CurSize.columns>(cur);
        } else if constexpr (sizeof...(LeftSizes) > 0){
            select_size_and_start<PPType, VType, VFlowType, LeftSizes...>(cur);
        } else{
            start_on_field<PPType, VType, VFlowType>(cur);
        }
    }

    template <class PPType, class VType, class VFlowType, sz... LeftSizes>
    void select_size_and_start_with_sizes(const Durex& cur) {
        if constexpr (sizeof...(LeftSizes) > 0) {
            select_size_and_start<PPType, VType, VFlowType, LeftSizes...>(cur);
        } else {
            start_on_field<PPType, VType, VFlowType>(cur);
        }
    }

    template <class CompTypesList, class OurTypesList, sz... AllSizes>
    class TypeSelector;

    template<class... CompTypes, class... OurTypes, sz... AllSizes>
    class TypeSelector<List<CompTypes...>, List<OurTypes...>, AllSizes...> {
        static constexpr bool double_flag = std::disjunction_v<std::is_same<CompTypes, DOUBLE>...>;
        static constexpr bool float_flag = std::disjunction_v<std::is_same<CompTypes, FLOAT>...>;

    public:
        template <class... DirtyTypes>
         static void select_type_and_size_and_start(const Durex& cur, std::string_view cur_type_name, DirtyTypes...types) {
                if(float_flag &&cur_type_name == kFloatTypeName) {
                    get_next_type_and_start<FLOAT, DirtyTypes...>(cur, types...);
                    return;
                }
                if(double_flag && cur_type_name == kDoubleTypeName) {
                    get_next_type_and_start<DOUBLE, DirtyTypes...>(cur, types...);
                    return;
                }

                std::size_t n{};
                std::size_t k{};
                if(parse_fixed(cur_type_name, FFpref, FFsuf, n, k)) {
                    if(approve_fixed_base<true, DirtyTypes...>(n, k, cur, types...)) {
                        return;
                    }
                }

                if(parse_fixed(cur_type_name, Fpref, Fsuf, n, k)) {
                    if(approve_fixed_base<false, DirtyTypes...>(n, k, cur, types...)) {
                        return;
                    }
                }

                throw std::invalid_argument("Invalid type name");
         }
    private:
        static bool parse_fixed(std::string_view cur_type_name, std::string_view pref,
                                std::string_view suf, std::size_t& n, std::size_t& k) {
            if(!cur_type_name.starts_with(pref) || !cur_type_name.ends_with(suf)) {
                return false;
            }
            cur_type_name.remove_prefix(pref.size());
            cur_type_name.remove_suffix(suf.size());

            const std::size_t pos = cur_type_name.find(',');

            if(pos == std::string_view::npos) {
                return false;
            }

            auto trim = [](std::string_view cur) {
                while(cur.starts_with(' ')) {
                    cur.remove_prefix(1);
                }
                while(cur.ends_with(' ')) {
                    cur.remove_suffix(1);
                }
                return cur;
            };

            const std::string_view new_n = trim(cur_type_name.substr(0, pos + 1));
            const std::string_view new_k = trim(cur_type_name.substr(pos + 1));

            if(std::from_chars(new_n.data(), new_n.data() + new_n.size(), n).ec != std::errc()) {
                return false;
            }

            if(std::from_chars(new_k.data(), new_k.data() + new_k.size(), k).ec != std::errc()) {
                return false;
            }

            return n > 0 && k > 0;
        }

        template <bool Fast, class... DirtyTypes>
        static bool approve_fixed_base(std::size_t n, std::size_t k, const Durex& cur, DirtyTypes... types) {
            return approve_fixed<0, Fast, DirtyTypes...>(n, k, cur, types...);
        }

        template <std::size_t ind, bool Fast, class... DirtyTypes>
        static bool approve_fixed(std::size_t n, std::size_t k, const Durex& cur, DirtyTypes...types) {
            using CurType = Get_by_ind_t<ind, CompTypes...>;

            if constexpr (requires {
                {
                    CurType::NFromFixed == std::size_t{}
                } -> std::same_as<bool>;

                {
                    CurType::KFromFixed == std::size_t{}
                } -> std::same_as<bool>;

                {
                    CurType::FastFromFixed == bool{}
                } -> std::same_as<bool>;
                }) {
                if constexpr (CurType::FastFromFixed == Fast) {
                    if(CurType::NFromFixed == n && CurType::KFromFixed == k) {
                        get_next_type_and_start<CurType, DirtyTypes...>(cur, types...);
                        return true;
                    }
                }
            }

            if constexpr (ind + 1 < sizeof...(CompTypes)) {
                return approve_fixed<ind + 1, Fast, DirtyTypes...>(n, k, cur, types...);
            } else {
                return false;
            }
        }

        template <class CurType, class... DirtyTypes>
        static void get_next_type_and_start(const Durex& cur, DirtyTypes... types) {
            if constexpr(sizeof...(types) > 0) {
                TypeSelector<List<CompTypes...>, List<OurTypes..., CurType>, AllSizes...>::
                select_type_and_size_and_start(cur, types...);
            } else {
                select_size_and_start_with_sizes<OurTypes..., CurType, AllSizes...>(cur);
            }
        }
    };

    struct RunTimeSettings {
        std::string_view p_type_name{};
        std::string_view v_type_name{};
        std::string_view v_flow_type_name{};
    };

    template <class List, sz... Sizes>
    class Go;

    template <class... FloatTypes, sz... AllSizes>
    class Go<List<FloatTypes...>, AllSizes...> {
    public:
        static Go SaveInfo(RunTimeSettings Settings) {
            return Go{std::move(Settings)};
        }

        void start_on_field(const Durex& cur) {
            TypeSelector<List<FloatTypes...>, List<>, AllSizes...>::
            select_type_and_size_and_start(cur, params.p_type_name, params.v_type_name, params.v_flow_type_name);
        }

    private:
        explicit constexpr Go(RunTimeSettings&& params) : params{std::move(params)} {}

        RunTimeSettings params;
    };
}

#ifndef TYPES
#error "TYPES is not defined"
#endif

#ifdef SIZES

using Go = types::Go<types::List<TYPES>, SIZES>;

#else

using Go = types::Go<types::List<TYPES>>;

#endif
void read_from_file_and_start(int argc = 0, char* argv[] = nullptr) {
    types::RunTimeSettings settings;

    std::string p_type, v_type, v_flow_type;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.starts_with("--p-type=")) {
            p_type = arg.substr(9);
            settings.p_type_name = p_type;
        } else if (arg.starts_with("--v-type=")) {
            v_type = arg.substr(9);
            settings.v_type_name = v_type;
        } else if (arg.starts_with("--v-flow-type=")) {
            v_flow_type = arg.substr(14);
            settings.v_flow_type_name = v_flow_type;
        }
    }

    Go simulator = Go::SaveInfo(settings);

    std::ifstream in("/mnt/d/Rcmple/cpp_project_2/example_field.txt");
    types::Durex cur;
    in >> cur.n >> cur.m;
    cur.resize();
    in >> cur.water_rho >> cur.air_rho >> cur.g;
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    for(auto& str : cur.field) {
        std::string cur_str;
        getline(in, cur_str);
        if(!cur_str.empty()) {
            if(cur_str.back() == '\r') {
                cur_str.pop_back();
            }
        }
        std::ranges::copy(cur_str, str.begin());
    }

    simulator.start_on_field(cur);
}
int main(int argc, char* argv[]) {
    bool p_type_set = false;
    bool v_type_set = false;
    bool v_flow_type_set = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.starts_with("--p-type=")) {
            p_type_set = true;
        } else if (arg.starts_with("--v-type=")) {
            v_type_set = true;
        } else if (arg.starts_with("--v-flow-type=")) {
            v_flow_type_set = true;
        }
    }

    if (p_type_set && v_type_set && v_flow_type_set) {
        read_from_file_and_start(argc, argv);
    } else {
        std::cerr << "Error: Missing required parameters. Please provide --p-type, --v-type, and --v-flow-type." << std::endl;
        return 1;
    }

    return 0;
}
