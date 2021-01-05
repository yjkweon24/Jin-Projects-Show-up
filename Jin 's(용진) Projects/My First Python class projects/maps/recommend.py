"""A Yelp-powered Restaurant Recommendation Program"""

from abstractions import *
from data import ALL_RESTAURANTS, CATEGORIES, USER_FILES, load_user_file
from ucb import main, trace, interact
from utils import distance, mean, zip, enumerate, sample
from visualize import draw_map

##################################
# Phase 2: Unsupervised Learning #
##################################


def find_closest(location, centroids):
    """Return the centroid in centroids that is closest to location.
    If multiple centroids are equally close, return the first one.

    >>> find_closest([3.0, 4.0], [[0.0, 0.0], [2.0, 3.0], [4.0, 3.0], [5.0, 5.0]])
    [2.0, 3.0]
    >>> find_closest([6, 1],[[1, 5], [3, 3]])
    [3,3]
    """
    # BEGIN Question 3
    "*** REPLACE THIS LINE ***"
    # For lists, dont need to use brace centroids[i] like this.
    return min(centroids, key = lambda i: distance(location, i))
    # END Question 3

def group_by_first(pairs):
    """Return a list of pairs that relates each unique key in the [key, value]
    pairs to a list of all values that appear paired with that key.

    Arguments:
    pairs -- a sequence of pairs

    >>> example = [ [1, 2], [3, 2], [2, 4], [1, 3], [3, 1], [1, 2] ]
    >>> group_by_first(example)
    [[2, 3, 2], [2, 1], [4]]
    """
    keys = []
    for key, _ in pairs:
        if key not in keys:
            keys.append(key)
    return [[y for x, y in pairs if x == key] for key in keys]


def group_by_centroid(restaurants, centroids):
    """Return a list of clusters, where each cluster contains all restaurants
    nearest to a corresponding centroid in centroids. Each item in
    restaurants should appear once in the result, along with the other
    restaurants closest to the same centroid.

    Q: If centroids is [[-1, 1], [5, -1], [1, 10], [-1, -10]],
    to which centroid will [6, 0] be associated?
    Choose the number of the correct choice:
    0) [-1, -10]
    1) [5, -1]
    2) [-1, 1]
    3) [1, 10]
    ? 1
    -- OK! --
    """
    # BEGIN Question 4
    "*** REPLACE THIS LINE ***"
    # We need to get location clusters first, and compare with 
    # centroids, and get the closest locations, and use group_by_first
    # to group by distances
    locations = [restaurant_location(i) for i in restaurants]
    closers = [find_closest(j, centroids) for j in locations]
    return group_by_first([[i, j] for i, j in zip(closers, restaurants)])
    # END Question 4


def find_centroid(cluster):
    """Return the centroid of the locations of the restaurants in cluster."""
    '''
    >>> cluster1 = [
    ...     make_restaurant('A', [-3, -4], [], 3, [make_review('A', 2)]),
    ...     make_restaurant('B', [1, -1],  [], 1, [make_review('B', 1)]),
    ...     make_restaurant('C', [2, -4],  [], 1, [make_review('C', 5)]),
    ... ]
    >>> find_centroid(cluster1) # should be a pair of decimals
    [0.0,-3.0]
    '''
    # BEGIN Question 5
    "*** REPLACE THIS LINE ***"
    locations_lat, locations_lon = [], []
    for i in cluster:
        #locations_lat += restaurant_location(i)[0]
        #locations_lon += restaurant_location(i)[1]
        locations_lat.append(restaurant_location(i)[0])
        locations_lon.append(restaurant_location(i)[1])
    return [mean(locations_lat), mean(locations_lon)]
    # END Question 5


def k_means(restaurants, k, max_updates=100):
    """Use k-means to group restaurants by location into k clusters."""
    '''
    Q: What is the first step of the iterative portion of the
    k-means algorithm?
    Choose the number of the correct choice:
    0) randomly initialize k centroids
    1) find the centroid (average position) of each cluster.
    2) create a cluster for each centroid consisting of all elements closest to
       that centroid.
    2
    '''
    assert len(restaurants) >= k, 'Not enough restaurants to cluster'
    old_centroids, n = [], 0
    # Select initial centroids randomly by choosing k different restaurants
    centroids = [restaurant_location(r) for r in sample(restaurants, k)]

    while old_centroids != centroids and n < max_updates:
        old_centroids = centroids
        # BEGIN Question 6
        "*** REPLACE THIS LINE ***"
        group_together = group_by_centroid(restaurants, old_centroids)
        centroids = [find_centroid(i) for i in group_together]
        # END Question 6
        n += 1
    return centroids


################################
# Phase 3: Supervised Learning #
################################


def find_predictor(user, restaurants, feature_fn):
    """Return a rating predictor (a function from restaurants to ratings),
    for a user by performing least-squares linear regression using feature_fn
    on the items in restaurants. Also, return the R^2 value of this model.

    Arguments:
    user -- A user
    restaurants -- A sequence of restaurants
    feature_fn -- A function that takes a restaurant and returns a number
    """
    '''
    Q: What does the list xs represent?
    Choose the number of the correct choice:
    0) the restaurants reviewed by user
    1) the restaurants in restaurants
    2) the extracted values for each restaurant in restaurants
    3) the names of restaurants in restaurants
    2
    Q: What does the list ys represent?
    Choose the number of the correct choice:
    0) the names for the restaurants in restaurants
    1) the average rating for the restaurants in restaurants
    2) the names for the restaurants reviewed by user
    3) user's ratings for the restaurants in restaurants
    3
    Q: Given a user, a list of restaurants, and a feature function, what
    does find_predictor from Problem 7 return?
    Choose the number of the correct choice:
    0) a predictor function, and its r_squared value
    1) a predictor function
    2) a restaurant
    3) an r_squared value
    0
    Q: After getting a list of [predictor, r_squared] pairs,
    which predictor should we select?
    Choose the number of the correct choice:
    0) an arbitrary predictor
    1) the predictor with the lowest r_squared
    2) the first predictor in the list
    3) the predictor with the highest r_squared
    3
    '''
    reviews_by_user = {review_restaurant_name(review): review_rating(review)
                       for review in user_reviews(user).values()}

    xs = [feature_fn(r) for r in restaurants]
    ys = [reviews_by_user[restaurant_name(r)] for r in restaurants]

    # BEGIN Question 7
    "*** REPLACE THIS LINE ***"
    # Initialization
    b, a, r_squared, S_xx, S_yy, S_xy = 0, 0, 0, 0, 0, 0
    combined = zip (xs, ys)
    mean_xs = mean(xs)
    mean_ys = mean(ys)
    for i in range(len(restaurants)):
        S_xx += (combined[i][0] - mean_xs) ** 2
        S_yy += (combined[i][1] - mean_ys) ** 2
        S_xy += (combined[i][0] - mean_xs) * (combined[i][1] - mean_ys)

    b = S_xy / S_xx
    a = mean(ys) - (b * mean(xs))
    r_squared = ((S_xy) ** 2) / (S_xx * S_yy)
    # END Question 7

    def predictor(restaurant):
        return b * feature_fn(restaurant) + a

    return predictor, r_squared

def best_predictor(user, restaurants, feature_fns):
    """Find the feature within feature_fns that gives the highest R^2 value
    for predicting ratings by the user; return a predictor using that feature.

    Arguments:
    user -- A user
    restaurants -- A list of restaurants
    feature_fns -- A sequence of functions that each takes a restaurant
    """
    '''
    Q: In best_predictor, what does the variable reviewed represent?
    Choose the number of the correct choice:
    0) a list of restaurants reviewed by the user
    1) a list of all possible restaurants
    2) a list of ratings for restaurants reviewed by the user
    0
    '''
    reviewed = user_reviewed_restaurants(user, restaurants)
    # BEGIN Question 8
    "*** REPLACE THIS LINE ***"
    max_feature = max(feature_fns, key = lambda x: find_predictor(user, reviewed, x)[1])
    return find_predictor(user, reviewed, max_feature)[0]
    # END Question 8


def rate_all(user, restaurants, feature_fns):
    """Return the predicted ratings of restaurants by user using the best
    predictor based on a function from feature_fns.

    Arguments:
    user -- A user
    restaurants -- A list of restaurants
    feature_fns -- A sequence of feature functions
    """
    '''
    Q: rate_all returns a dictionary. What are the keys of this dictionary?
    Choose the number of the correct choice:
    0) restaurants
    1) restaurant ratings
    2) restaurant names
    2
    Q: What are the values of the returned dictionary?
    Choose the number of the correct choice:
    0) numbers - user ratings only
    1) numbers - mean restaurant ratings
    2) lists - list of all restaurant ratings
    3) numbers - predicted ratings only
    4) numbers - a mix of user ratings and predicted ratings
    4
    Q: In rate_all, what does the variable reviewed represent?
    Choose the number of the correct choice:
    0) a list of ratings for restaurants reviewed by the user
    1) a list of all possible restaurants
    2) a list of restaurants reviewed by the user
    2

    '''
    predictor = best_predictor(user, ALL_RESTAURANTS, feature_fns)
    reviewed = user_reviewed_restaurants(user, restaurants)
    # BEGIN Question 9
    "*** REPLACE THIS LINE ***"
    name_and_ratings = {}
    for i in restaurants:
        if i in reviewed:
            name_and_ratings[restaurant_name(i)] = user_rating(user, restaurant_name(i))
        else:
            name_and_ratings[restaurant_name(i)] = predictor(i)
    return name_and_ratings
    # END Question 9

def search(query, restaurants):
    """Return each restaurant in restaurants that has query as a category.

    Arguments:
    query -- A string
    restaurants -- A sequence of restaurants
    """
    '''
    Q: When does a restaurant match a search query?
    Choose the number of the correct choice:
    0) if the query string is one of the restaurant's categories
    1) if the query string is a substring of the restaurant's name
    2) if the query string is equal to the restaurant's categories
    3) if the query string is mentioned in the restaurant's reviews
    0
    Q: What type of object does search return?
    0) a dictionary that maps restaurant categories (strings) to restaurants
    1) a list of restaurant names (strings)
    2) a dictionary that maps restaurant names (strings) to restaurants
    3) a list of restaurants
    3
    '''
    # BEGIN Question 10
    "*** REPLACE THIS LINE ***"
    return [x for x in restaurants if query in restaurant_categories(x)] 
    # END Question 10

def feature_set():
    """Return a sequence of feature functions."""
    return [lambda r: mean(restaurant_ratings(r)),
            restaurant_price,
            lambda r: len(restaurant_ratings(r)),
            lambda r: restaurant_location(r)[0],
            lambda r: restaurant_location(r)[1]]


@main
def main(*args):
    import argparse
    parser = argparse.ArgumentParser(
        description='Run Recommendations',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-u', '--user', type=str, choices=USER_FILES,
                        default='test_user',
                        metavar='USER',
                        help='user file, e.g.\n' +
                        '{{{}}}'.format(','.join(sample(USER_FILES, 3))))
    parser.add_argument('-k', '--k', type=int, help='for k-means')
    parser.add_argument('-q', '--query', choices=CATEGORIES,
                        metavar='QUERY',
                        help='search for restaurants by category e.g.\n'
                        '{{{}}}'.format(','.join(sample(CATEGORIES, 3))))
    parser.add_argument('-p', '--predict', action='store_true',
                        help='predict ratings for all restaurants')
    parser.add_argument('-r', '--restaurants', action='store_true',
                        help='outputs a list of restaurant names')
    args = parser.parse_args()

    # Output a list of restaurant names
    if args.restaurants:
        print('Restaurant names:')
        for restaurant in sorted(ALL_RESTAURANTS, key=restaurant_name):
            print(repr(restaurant_name(restaurant)))
        exit(0)

    # Select restaurants using a category query
    if args.query:
        restaurants = search(args.query, ALL_RESTAURANTS)
    else:
        restaurants = ALL_RESTAURANTS

    # Load a user
    assert args.user, 'A --user is required to draw a map'
    user = load_user_file('{}.dat'.format(args.user))

    # Collect ratings
    if args.predict:
        ratings = rate_all(user, restaurants, feature_set())
    else:
        restaurants = user_reviewed_restaurants(user, restaurants)
        names = [restaurant_name(r) for r in restaurants]
        ratings = {name: user_rating(user, name) for name in names}

    # Draw the visualization
    if args.k:
        centroids = k_means(restaurants, min(args.k, len(restaurants)))
    else:
        centroids = [restaurant_location(r) for r in restaurants]
    draw_map(centroids, restaurants, ratings)