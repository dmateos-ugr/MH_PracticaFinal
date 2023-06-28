const std = @import("std");
const Random = std.rand.Random;
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;

const SortContext = struct {
    population: [][]const f64,
    fitnesses: []f64,

    pub fn lessThan(self: SortContext, idx1: usize, idx2: usize) bool {
        return self.fitnesses[idx1] < self.fitnesses[idx2];
    }

    pub fn swap(self: *SortContext, idx1: usize, idx2: usize) void {
        std.mem.swap([]const f64, &self.population[idx1], &self.population[idx2]);
        std.mem.swap(f64, &self.fitnesses[idx1], &self.fitnesses[idx2]);
    }
};

pub fn sortPopulationBestFirst(population: [][]const f64, fitnesses: []f64) void {
    var sort_context = SortContext{ .population = population, .fitnesses = fitnesses };
    std.sort.insertionSortContext(population.len, &sort_context);
    std.debug.assert(std.sort.isSorted(f64, fitnesses, {}, std.sort.asc(f64)));
}


pub fn randomFloat(min: f64, max: f64, rnd: Random) f64 {
    assert(min <= max);
    return min + rnd.float(f64)*(max - min);
}

pub fn createRandomSolution(n: usize, allocator: Allocator, rnd: Random) ![]f64 {
    const sol = try allocator.alloc(f64, n);
    for (sol) |*value| {
        value.* = randomFloat(-100, 100, rnd);
    }
    return sol;
}

pub fn euclideanNorm(v: []f64) f64 {
	var length: f64 = 0;
	for (v) |vi| {
		length += vi * vi;
	}
	return std.math.sqrt(length);
}

// pub fn vecAdd

pub fn distance(v1: []const f64, v2: []const f64) f64 {
	var dist: f64 = 0;
	for (v1, v2) |e1, e2| {
		const diff = e1 - e2;
		dist += diff * diff;
	}
	return std.math.sqrt(dist);
}

pub fn normalize(v: []f64) void {
	const length = euclideanNorm(v);
	if (length == 0)
		return;
	for (v) |*vi| {
		vi.* /= length;
		std.debug.assert(!std.math.isInf(vi.*));
	}
}

pub fn randomDirection(direction: []f64, rnd: Random) void {
	for (direction) |*dir_i| {
		dir_i.* = rnd.float(f64);
	}
	normalize(direction);
}
