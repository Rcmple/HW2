
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
            EType &add(int x, int y, int dx, int dy, EType& dv) {
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

        explicit constexpr Simulator(const Durex& cur) noexcept requires (static_flag) {
            for(std::size_t ind = 0; auto& str : field) {
                const std::vector<char>& cur_str = cur.field[ind];
                assert(cur_str.size() == str.size());
                std::ranges::copy(cur_str, str.begin());
            }
        }
        explicit constexpr Simulator(const Durex& cur)
        requires(!static_flag)
                : Simulator(cur, cur.field.size(), cur.field.front().size() - 1) {}
        [[nodiscard]] constexpr std::size_t get_n() const noexcept {
            return field.size();
        }
        [[nodiscard]] constexpr std::size_t get_m() const noexcept {
            return field.front().size() - 1;
        }

        void start() {

        }
    private:
        explicit constexpr Simulator(const Durex& cur, std::size_t rows, std::size_t columns) requires(!static_flag):
                field{cur.field},
                p{rows, typename PPField::value_type(columns)},
                p_old{rows, typename PPField::value_type(columns)},
                v_field{rows, typename VField::CurType(columns)},
                v_flow_field{rows, typename VFlowField::CurType(columns)},
                last_use_field{rows, typename LastUseField::value_type(columns)} {}
        Field field{};
        PPField p, p_old{};
        VField v_field{};
        VFlowField v_flow_field{};
        LastUseField last_use_field{};
        int UT = 0;
        std::array <PPType, 256> rho{};
    };

    template <class PPType, class VType, class VFlowType,
            std::size_t Rows = dyn_flag, std::size_t Columns = dyn_flag>
    void start_on_field(const Durex& cur) {
        Simulator<PPType, VType, VFlowType, Rows, Columns> a{cur};
        a.start();
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
