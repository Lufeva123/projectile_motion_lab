import math
from decimal import Decimal
import sys

g = 9.807
DECIMAL_PLACES = 4


def calculate_delta_h (h1, h2):
    return round(h2 - h1, DECIMAL_PLACES)


def calculate_delta_h_prime(h1_prime, h2_prime):
    return round(h2_prime - h1_prime, DECIMAL_PLACES)


def calculate_v_0 (delta_h, delta_h_prime):
    return round(math.sqrt(10/7 * g * (delta_h_prime - delta_h)), DECIMAL_PLACES)


def calculate_v_0_x(v_0, d, l):
    return round((v_0 * d)/l, DECIMAL_PLACES)


def calculate_v_0_y(v_0, h2, h3, l):
    return round((v_0 * (h2 - h3))/l, DECIMAL_PLACES)


def calculate_t(v_0_y, h2):
    delta = (v_0_y ** 2) + (2 * g * h2)
    return round((v_0_y + math.sqrt(delta))/g, DECIMAL_PLACES)


def calculate_range(v_0_x, t):
    return round(v_0_x * t, DECIMAL_PLACES)


if __name__ == "__main__":
    h1 = float(input("h1 (meters) = "))
    h1_prime = float(input("h1_prime (meters) = "))
    h2 = float(input("h2 (meters) = "))
    h2_prime = float(input("h2_prime (meters) = "))
    h3 = float(input("h3 (meters) = "))
    d = float(input("d (meters) = "))
    l = float(input("l (meters) = "))

    delta_h = calculate_delta_h(h1, h2)
    delta_h_prime = calculate_delta_h_prime(h1_prime, h2_prime)

    v_0 = calculate_v_0(delta_h, delta_h_prime)

    v_0_x = calculate_v_0_x(v_0, d, l)
    v_0_y = calculate_v_0_y(v_0, h2, h3, l)

    t = calculate_t(v_0_y, h2)

    rng = calculate_range(v_0_x, t)

    print()
    print("Range = {:.5f} meters".format(rng))
    print()
