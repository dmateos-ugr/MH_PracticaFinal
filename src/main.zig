const std = @import("std");
const utils = @import("utils.zig");
const c = @cImport({
    @cInclude("cec2017/cec17.h");
});
const Allocator = std.mem.Allocator;
const Random = std.rand.Random;
const assert = std.debug.assert;

const N_POPULATION = 20;

fn ppa(funcid: usize, n: usize, allocator: Allocator, rnd: Random) ![]f64 {
    c.cec17_init("ppa", @intCast(c_int, funcid), @intCast(c_int, n));

    const max_evaluations = n * 10000;
    var population: [N_POPULATION][]f64 = undefined;
    var fitnesses: [N_POPULATION]f64 = undefined;

    // Allocate memory for the population, initializing it randomly
    for (&population) |*element| {
        element.* = try utils.createRandomSolution(n, allocator, rnd);
    }
    defer for (population) |element| {
        allocator.free(element);
    };

    var evaluations: usize = 0;
    while (evaluations < max_evaluations) {
        // Get fitnesses, and sort population according to them
        for (population, &fitnesses) |element, *fitness| {
            fitness.* = c.cec17_fitness(element.ptr);
        }
        utils.sortPopulationBestFirst(&population, &fitnesses);

    }

    return allocator.dupe(f64, population[0]);
}

const PREDATOR_IDX = N_POPULATION - 1;
const BEST_PREY_IDX = 0;
const K = 50;
const LAMBDA_MAX = 7;
const LAMBDA_MIN = 3;
const PROB_FOLLOW_UP = 0.5;
const TAU = 0.3;
fn movePrey(
    idx: usize,
    population: [][]f64,
    fitnesses: []f64,
    allocator: Allocator,
    rnd: Random,
) usize {
    const n = fitnesses.len;

    // fitnesses ordenado de mayor a menor
    const num_preys_with_better_fitness = for (fitnesses, 0..) |fitness, i| {
        if (fitness <= fitnesses[idx])
            break i;
    } else 0;

    var direction = try allocator.alloc(f64, n);
    defer allocator.free(direction);
    var mut = try allocator.alloc(f64, n);
    defer allocator.free(mut);

    if (num_preys_with_better_fitness == 0) {
        // best prey (puede haber varias)
        // generate K random directions
        // perform local search with K mutations
        for (0..K) |_| {
            // Calculate mutation
            utils.randomDirection(direction, rnd);
            for (0..n) |i| {
                mut[i] = population[idx][i] + LAMBDA_MIN * rnd.float(f64) * direction[i];
            }

            const mut_fitness = c.cec17_fitness(mut.ptr);
            if (mut_fitness > fitnesses[idx])
                @memcpy(population[idx], mut);
        }
        return K;
    }

    if (rnd.float() <= PROB_FOLLOW_UP) {
        // TODO este doble bucle esta diferente a pabloco
        for (0..n) |i| {
            direction[i] = 0;
            for (0..num_preys_with_better_fitness) |jdx| {
                const dist = utils.distance(population[idx], population[jdx]);
                const tmp = std.math.pow(f64, fitnesses[jdx], TAU);
                direction[i] += std.math.exp(tmp - dist) * (population[jdx][i] - population[idx][i]);
            }
        }
        // TODO terminar esta parte

    } else {
        utils.randomDirection(direction, rnd);

        for (0..n) |i| {
            mut[i] = population[idx][i] + direction[i];
        }
        const d1 = utils.distance(population[PREDATOR_IDX], mut);

        for (0..n) |i| {
            mut[i] = population[idx][i] - direction[i];
        }
        const d2 = utils.distance(population[PREDATOR_IDX], mut);

        const mult = if (d1 <= d2) -1 else 1;
        for (0..n) |i| {
            population[idx][i] += mult * LAMBDA_MAX * rnd.float(f64) * direction[i];
        }
        return 0;
    }
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const seed = 1234;
    var rng = std.rand.DefaultPrng.init(seed);
    const rnd = rng.random();

    for (1..31) |funcid| {
        const sol = try ppa(funcid, 10, allocator, rnd);
        defer allocator.free(sol);
    }
}
