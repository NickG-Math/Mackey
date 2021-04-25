#include "../BC4S2.hpp"

namespace mackey
{
	template<typename group_t>
	BC4S2<group_t>::BC4S2(int dimension) : dimension(dimension)
	{
		cells_equivariant.resize(dimension + 1);
		cuts.resize(dimension + 1);
		cell_compressed init;
		init.cell.resize(dimension + 1);
		std::iota(init.cell.begin(), init.cell.end(), 0);
		cells_equivariant[dimension] = { init };
		for (int i = 0; i < dimension; i++)
		{
			auto pair = quickcut(cells_equivariant[dimension - i]);
			cells_equivariant[dimension - i - 1] = pair.first;
			cuts[dimension - i] = pair.second;
		}
		setcells();
		setexpandedcells();
	}

	template<typename group_t>
	int BC4S2<group_t>::getrank(const cell_compressed& v) const {
		for (const auto i : v.cell)
			if (i % 4 == 2 || i % 4 == 3) {
				//check if it's only 2,3
				for (const auto j : v.cell)
					if (j % 4 == 0 || j % 4 == 1)
						return 4;
				return 2;
			}
		char oddoreven = v.cell[0] % 4;
		for (const auto i : v.cell) {
			if (i % 4 != oddoreven)
				return 2;
		}
		return 1;
	}

	template<typename group_t>
	auto BC4S2<group_t>::getrank(const std::vector<cell_compressed>& v) const -> rank_t {
		rank_t r(v.size());
		int counter = 0;
		for (size_t i = 0; i < v.size(); i++)
			r[i] = getrank(v[i]);
		return r;
	}

	template<typename group_t>
	void BC4S2<group_t>::setcells() {
		this->cells.resize(dimension + 1);
		for (int i = 0; i <= dimension; i++)
			this->cells[i] = getrank(cells_equivariant[i]);
	}

	template<typename group_t>
	void BC4S2<group_t>::setexpandedcells() {
		cells_nonequivariant.resize(dimension + 1);
		for (int i = 0; i <= dimension; i++)
			for (int j = 0; j < cells_equivariant[i].size(); j++) {
				cells_nonequivariant[i].push_back(getexpandedcell(i, eq_t(j)));
				for (auto k = 1; k < this->cells[i][j]; k++)
					cells_nonequivariant[i].push_back(cells_nonequivariant[i].back().action());
			}
	}

	template<typename group_t>
	auto BC4S2<group_t>::getfaces(int dim, eq_t k) -> std::vector<cell_expanded> {
		std::vector<cell_expanded> faces;
		const auto totalcuts = cuts[dim][k.index].size();
		faces.reserve(2 * totalcuts);
		auto cell = cells_nonequivariant[dim][this->convert(k, dim).index];
		for (size_t i = 0; i < totalcuts; i++) {
			auto cut = cuts[dim][k.index][i];
			auto face = cell;
			face.cell[cut] = 0;
			faces.push_back(face);
			if (cut % 4 == 2 && (i + 1 == totalcuts || cuts[dim][k.index][i + 1] != cut + 1) && cell.cell[cut + 1] == 1) {
				auto face = cell;
				face.cell[cut + 1] = 0;
				faces.push_back(face);
			}
			if (cut % 4 == 3 && (i == 0 || cuts[dim][k.index][i - 1] != cut - 1) && cell.cell[cut - 1] == 1) {
				auto face = cell;
				face.cell[cut - 1] = 0;
				faces.push_back(face);
			}
		}
		return faces;
	}

	template<typename group_t>
	auto BC4S2<group_t>::getexpandedcell(int dim, eq_t i) const -> cell_expanded {
		cell_expanded ce;
		ce.cell.resize(dimension + 1);
		if (dim == 0) {
			ce.cell[cells_equivariant[dim][i.index].cell[0]] = 1;
			return ce;
		}
		for (const auto j : cuts[dim][i.index]) {
			ce.cell[j] = 1;
			if (j % 4 == 2)
				ce.cell[j + 1] = 3; //tentative, have to check if cells_equivariant contain j+1
			if (j % 4 == 3)
				ce.cell[j - 1] = 3; //tentative, have to check if cells_equivariant contain j+1
		}
		for (const auto j : cells_equivariant[dim][i.index].cell) {
			if (ce.cell[j] == 0)
				ce.cell[j] = 2; //placeholder
			if (ce.cell[j] == 3)
				ce.cell[j] = 1; //fix the 3
		}
		for (auto& j : ce.cell)
			if (j == 3) //fix the 3 in the case where it's not found
				j = 0;
		return ce;
	}

	template<typename group_t>
	int64_t BC4S2<group_t>::estimate_nonzero_entries(int i)
	{
		return -1;
	}

	template<typename group_t>
	auto BC4S2<group_t>::boundary(int dim, eq_t k)
	{
		static_assert(scalar_t<diff_t>::order == 2, "Currently only supported for Z/2 coefficients");

		boundary_t bd;
		std::vector<scalar_t<diff_t>> coeff_per_cell(cells_nonequivariant[dim - 1].size());
		auto faces = getfaces(dim, k);
		for (const auto& face : faces) {

			//a face with a 2 may need to split to +1 and -1 if the subspaces have no 2 there
			//we record where the splits must happen
			std::vector<int> splits;
			for (size_t j = 0; j < face.cell.size(); j++) {
				if (face.cell[j] == 2) //check if any cell here uses 1 instead of 2, then split. 
					for (const auto& lowcell : cells_nonequivariant[dim - 1]) {
						if (lowcell.isfaceof(face) && lowcell.cell[j] == 1) {
							splits.push_back(j);
							break;
						}
					}
			}

			//now we split based on all the combinations
			for (int j = 0; j <= splits.size(); j++) {
				auto cg = CombinationGenerator<int>(splits.size(), j);
				for (const auto& i : cg) {
					auto splitface = face;
					for (auto k : i)
						splitface.cell[splits[k]] = 1;
					for (auto& k : splitface.cell)
						if (k == 2)
							k = -1;
					//now we add the splitface to what we already have
					int ind = std::find(cells_nonequivariant[dim - 1].begin(), cells_nonequivariant[dim - 1].end(), splitface) - cells_nonequivariant[dim - 1].begin();
					coeff_per_cell[ind] += 1;
				}
			}
		}
		for (size_t i = 0; i < coeff_per_cell.size(); i++)
			if (coeff_per_cell[i] != 0)
				bd.push_back({ neq_t(i), coeff_per_cell[i] });
		return bd;
	}

	template<typename group_t>
	auto BC4S2<group_t>::quickcut(const cell_compressed& v) const ->std::array<std::pair<cell_compressed, int>, 2>{
		std::pair<cell_compressed, int> first, second;
		first = erase_by_criterion(v, [](int i)->bool {return i % 4 == 2 || i % 4 == 3;});
		if (!first.first.cell.empty())
			second = erase_by_criterion(v, [](int i)->bool {return i % 4 == 0 || i % 4 == 1;});
		else {
			first = erase_by_criterion(v, [](int i)->bool {return i % 4 == 0 || i % 4 == 1;});
			if (!first.first.cell.empty())
				second = erase_by_criterion(v, [=](int i)->bool {return i % 4 != first.second % 4;});
		}
		return { first,second };
	}

	template<typename group_t>
	auto BC4S2<group_t>::quickcut(const std::vector<cell_compressed>& v) const -> std::pair<std::vector<cell_compressed>, std::vector<std::vector<int>>> {
		std::vector<cell_compressed> newcells;
		std::vector<std::vector<int>> cut(v.size());
		for (auto s = 0; s < v.size(); s++) {
			auto pair = quickcut(v[s]);
			for (int t = 0; t <= 1; t++) {
				if (!pair[t].first.cell.empty() && std::find(newcells.begin(), newcells.end(), pair[t].first) == newcells.end())
					newcells.push_back(pair[t].first);
				if (!pair[t].first.cell.empty())
					cut[s].push_back(pair[t].second);
			}
		}
		return std::make_pair(newcells, cut);
	}

	///<Vector of 0,1 entries representing the homogeneous coordinates
	template<typename group_t>
	auto BC4S2<group_t>::cell_expanded::action() const -> cell_expanded {
		auto gx = cell;
		for (size_t i = 0; i < cell.size(); i++)
		{
			if (i % 4 == 1)
				gx[i] = -cell[i]; //antipode
			else if (i % 4 == 2)
			{ //rotation
				gx[i] = -cell[i + 1];
				gx[i + 1] = cell[i];
			}
		}
		return cell_expanded{ gx };
	}


	template<typename group_t>
	bool BC4S2<group_t>::cell_expanded::operator==(const cell_expanded& c) const {
		bool equal = 1;
		for (size_t i = 0; i < cell.size(); i++)
			if (abs(cell[i]) == 1 && abs(c.cell[i]) == 1) {//find first index non 0,2 to see whether cell=c or cell=-c
				equal = (cell[i] == c.cell[i]);
				break;
			}
		for (size_t i = 0; i < cell.size(); i++) {
			//if exactly one is 0 then bad
			if (cell[i] == 0 ^ c.cell[i] == 0)
				return 0;
			//either both are 0 or both nonzero
			//if both 0 then good
			else if (cell[i] == 0 && c.cell[i] == 0)
				continue;
			//both nonzero
			//if one is 2 then good
			else if (abs(cell[i]) == 2 || abs(c.cell[i]) == 2)
				continue;
			//both are pm 1
			else if ((equal && cell[i] != c.cell[i]) || (!equal && cell[i] != -c.cell[i]))
				return 0;
			//equal or -equal
		}
		return 1;
	}

	template<typename group_t>
	bool BC4S2<group_t>::cell_expanded::isfaceof(const cell_expanded& c) const {
		for (size_t i = 0; i < cell.size(); i++)
			if (cell[i] != 0 && c.cell[i] == 0)
				return 0;
		return 1;
	}


	template<typename group_t>
	bool BC4S2<group_t>::cell_compressed::operator==(const cell_compressed& c) const{
		return cell == c.cell;
	}

	template<typename group_t>
	template<typename T>
	auto BC4S2<group_t>::erase_by_criterion(const cell_compressed& v, const T& criterion) const -> std::pair<cell_compressed, int> {
		cell_compressed w;
		int val = -1;
		for (int i = v.cell.size() - 1; i >= 0; i--) {
			if (criterion(v.cell[i])) {
				w = v;
				w.cell.erase(w.cell.begin() + i);
				val = v.cell[i];
				break;
			}
		}
		return std::make_pair(w, val);
	}


}
