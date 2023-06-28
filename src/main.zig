const std = @import("std");
const utils = @import("utils.zig");
const c = @cImport({
    @cInclude("cec2017/cec17.h");
});
const Allocator = std.mem.Allocator;
const Random = std.rand.Random;
const assert = std.debug.assert;
const print = std.debug.print;

const N_POPULATION = 20;
const MEMETIC = false;
const N_EVALUATIONS_MEMETIC = 10;

fn surviveValue(solution: []f64) f64 {
    for (solution) |v| {
        std.debug.assert(!std.math.isInf(v));
        std.debug.assert(!std.math.isNan(v));
    }
    const fitness = c.cec17_fitness(solution.ptr);
    // return fitness;
    return std.math.log10(fitness);
    // return 1/std.math.log2(fitness);
}

fn busquedaLocal(
    evaluations: usize,
    element: []f64,
    fitness: *f64,
    allocator: Allocator,
    rnd: Random,
) !void {
    const n = element.len;
    const mut = try allocator.alloc(f64, n);
    defer allocator.free(mut);

    for (0..evaluations) |_| {
        utils.randomDirection(mut, rnd);
        const length = rnd.float(f64) * BL_MAX_STEP_LENGTH;
        for (0..n) |i| {
            mut[i] = element[i] + mut[i] * length;
        }
        const mut_fitness = surviveValue(mut);
        if (mut_fitness < fitness.*) {
            @memcpy(element, mut);
            fitness.* = mut_fitness;
        }
    }
}

fn ppa(funcid: usize, n: usize, allocator: Allocator, rnd: Random) ![]f64 {
    c.cec17_init("ppa", @intCast(c_int, funcid), @intCast(c_int, n));
    // c.cec17_print_output();

    const max_evaluations = n * 10000;
    var population1: [N_POPULATION][]f64 = undefined;
    var population2: [N_POPULATION][]f64 = undefined;
    var fitnesses1: [N_POPULATION]f64 = undefined;
    var fitnesses2: [N_POPULATION]f64 = undefined;

    // Allocate memory for the population, initializing it randomly
    for (&population1, &population2) |*element1, *element2| {
        element1.* = try utils.createRandomSolution(n, allocator, rnd);
        element2.* = try allocator.alloc(f64, n);
    }
    defer for (population1, population2) |element1, element2| {
        allocator.free(element1);
        allocator.free(element2);
    };

    var population = &population1;
    var fitnesses = &fitnesses1;
    var new_population = &population2;
    var new_fitnesses = &fitnesses2;

    // Get fitnesses, and sort population according to them
    for (population, fitnesses) |element, *fitness| {
        fitness.* = surviveValue(element);
    }
    // for (fitnesses) |fitness| {
    //     std.debug.print("{}\n", .{fitness});
    // }

    var evaluations: usize = 0;
    while (evaluations < max_evaluations) {
        // std.debug.print("evaluations: {}\n", .{evaluations});
        utils.sortPopulationBestFirst(population, fitnesses);

        // Move everything, updating new_population
        for (0..N_POPULATION - 1) |idx| {
            evaluations += try movePrey(
                idx,
                population,
                fitnesses,
                new_population[idx],
                &new_fitnesses[idx],
                allocator,
                rnd,
            );
        }

        evaluations += try movePredator(
            population,
            new_population[PREDATOR_IDX],
            &new_fitnesses[PREDATOR_IDX],
            allocator,
            rnd,
        );

        if (MEMETIC) {
            for (new_population, new_fitnesses) |element, *fitness| {
                try busquedaLocal(N_EVALUATIONS_MEMETIC, element, fitness, allocator, rnd);
                evaluations += N_EVALUATIONS_MEMETIC;
                if (evaluations >= max_evaluations)
                    break;
            }
        }

        std.mem.swap(*[N_POPULATION][]f64, &population, &new_population);
        std.mem.swap(*[N_POPULATION]f64, &fitnesses, &new_fitnesses);
    }

    return allocator.dupe(f64, population[0]);
}

const K = 20;
const PREDATOR_IDX = N_POPULATION - 1;
const BEST_PREY_IDX = 0;
const WORST_PREY_IDX = N_POPULATION - 2;
const BL_MAX_STEP_LENGTH = 15;
const LAMBDA_MAX = 7;
const LAMBDA_MIN = 3;
const PROB_FOLLOW_UP = 0.8;
const TAU = 0.8; //0.3;
const BETA = 1;


fn abs(v: f64) f64 {
    return if (v >= 0) v else -v;
}

fn movePrey(
    idx: usize,
    population: []const []const f64,
    fitnesses: []const f64,
    result: []f64,
    result_fitness: *f64,
    allocator: Allocator,
    rnd: Random,
) !usize {
    const n = population[0].len;
    @memcpy(result, population[idx]);
    result_fitness.* = fitnesses[idx];

    // fitnesses ordenado de menor a mayor
    std.debug.assert(std.sort.isSorted(f64, fitnesses, {}, std.sort.asc(f64)));
    const num_preys_with_better_fitness = for (fitnesses, 0..) |fitness, i| {
        if (fitness >= fitnesses[idx])
            break i;
    } else 0;

    const mut = try allocator.alloc(f64, n);
    defer allocator.free(mut);

    if (num_preys_with_better_fitness == 0) {
        // best prey (puede haber varias)
        try busquedaLocal(K, result, result_fitness, allocator, rnd);
        return K;
    }

    const direction = try allocator.alloc(f64, n);
    defer allocator.free(direction);
    if (rnd.float(f64) <= PROB_FOLLOW_UP) {
        for (0..n) |i| {
            direction[i] = 0;
            for (0..num_preys_with_better_fitness) |jdx| {
                const dist = utils.distance(population[idx], population[jdx]);
                const tmp = std.math.pow(f64, fitnesses[jdx], TAU);
                direction[i] += std.math.exp(tmp - dist) * (population[jdx][i] - population[idx][i]);
                assert(!std.math.isInf(direction[i]));
            }
        }

        const best_direction_bl = blk: {
            const direction_bl = try allocator.alloc(f64, n);
            defer allocator.free(direction_bl);
            const best_direction_bl = try allocator.alloc(f64, n);
            @memset(best_direction_bl, 0);
            var best_direction_bl_fitness = fitnesses[idx];

            for (0..K) |_| {
                // Calculate mutation
                utils.randomDirection(direction_bl, rnd);
                for (0..n) |i| {
                    mut[i] = population[idx][i] + LAMBDA_MIN * direction_bl[i];
                }

                const mut_fitness = surviveValue(mut);
                if (mut_fitness < best_direction_bl_fitness) {
                    @memcpy(best_direction_bl, direction_bl);
                    best_direction_bl_fitness = mut_fitness;
                }
            }
            break :blk best_direction_bl;
        };
        defer allocator.free(best_direction_bl);

        // Ya tenemos best_direction_bl (yr) y direction (yi). Calculamos el resultado.
        // print("{}\n", .{direction[0]});
        utils.normalize(direction);
        const lmax = LAMBDA_MAX / std.math.exp(BETA * abs(fitnesses[idx] - fitnesses[PREDATOR_IDX]));
        // const epsilon1 = rnd.float(f64);
        // const epsilon2 = rnd.float(f64);
        const epsilon1 = if (MEMETIC) 0.5 else 1;
        // const epsilon1 = 0.5;
        const epsilon2 = 1 - epsilon1;
        for (0..n) |i| {
            result[i] += lmax * epsilon1 * direction[i] + epsilon2 * best_direction_bl[i];
        }
        result_fitness.* = surviveValue(result);
        return K + 1;
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

        const mult: f64 = if (d1 <= d2) -1 else 1;
        for (0..n) |i| {
            result[i] += mult * LAMBDA_MAX * rnd.float(f64) * direction[i];
        }
        result_fitness.* = surviveValue(result);
        return 1;
    }
}

fn movePredator(
    population: []const []const f64,
    result: []f64,
    result_fitness: *f64,
    allocator: Allocator,
    rnd: Random,
) !usize {
    const n = population[0].len;

    const direction_rnd = try allocator.alloc(f64, n);
    defer allocator.free(direction_rnd);
    utils.randomDirection(direction_rnd, rnd);

    const direction_prey = try allocator.alloc(f64, n);
    defer allocator.free(direction_prey);
    for (0..n) |i| {
        direction_prey[i] = population[WORST_PREY_IDX][i] - population[PREDATOR_IDX][i];
    }
    utils.normalize(direction_prey);

    const epsilon1 = rnd.float(f64);
    const epsilon2 = rnd.float(f64);
    @memcpy(result, population[PREDATOR_IDX]);
    for (0..n) |i| {
        result[i] += LAMBDA_MAX * epsilon1 * direction_rnd[i] + LAMBDA_MIN * epsilon2 * direction_prey[i];
    }
    result_fitness.* = surviveValue(result);
    return 1;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const seed = 1234;
    var rng = std.rand.DefaultPrng.init(seed);
    const rnd = rng.random();

    for (1..31) |funcid| {
        print("funcid: {}\n", .{funcid});
        const sol = try ppa(funcid, 10, allocator, rnd);
        defer allocator.free(sol);
        print("error: {}\n", .{c.cec17_error(c.cec17_fitness(sol.ptr))});
    }
}
