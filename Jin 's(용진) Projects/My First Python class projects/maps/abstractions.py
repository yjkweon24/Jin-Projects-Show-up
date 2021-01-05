"""Data Abstractions"""

from utils import mean

#############################
# Phase 1: Data Abstraction #
#############################


# Reviews

def make_review(restaurant_name, rating):
    """Return a review data abstraction."""
    return [restaurant_name, rating]

def review_restaurant_name(review):
    """Return the restaurant name of the review, which is a string."""
    return review[0]

def review_rating(review):
    """Return the number of stars given by the review, which is a
    floating point number between 1 and 5."""
    return review[1]


# Users

def make_user(name, reviews):
    """Return a user data abstraction."""
    return [name, {review_restaurant_name(r): r for r in reviews}]

def user_name(user):
    """Return the name of the user, which is a string."""
    return user[0]

def user_reviews(user):
    """Return a dictionary from restaurant names to reviews by the user."""
    return user[1]


### === +++ USER ABSTRACTION BARRIER +++ === ###

def user_reviewed_restaurants(user, restaurants):
    """Return the subset of restaurants reviewed by user.

    Arguments:
    user -- a user
    restaurants -- a list of restaurant data abstractions
    """
    names = list(user_reviews(user))
    return [r for r in restaurants if restaurant_name(r) in names]

def user_rating(user, restaurant_name):
    """Return the rating given for restaurant_name by user."""
    reviewed_by_user = user_reviews(user)
    user_review = reviewed_by_user[restaurant_name]
    return review_rating(user_review)


# Restaurants

def make_restaurant(name, location, categories, price, reviews):
    """Return a restaurant data abstraction containing the name, location,
    categories, price, and reviews for that restaurant."""
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    # Make a dictionary
    return {"Name": name, "Location": location, "Categories": categories,
            "Price": price, "Reviews": [review_rating(x) for x in reviews]
        }
    # END Question 2

def restaurant_name(restaurant):
    """Return the name of the restaurant, which is a string."""
    '''
    >>> soda_reviews = [make_review('Soda', 4.5), make_review('Soda', 4)]
    >>> soda = make_restaurant('Soda', [127.0, 0.1],
                            ['Restaurants', 'Breakfast & Brunch'],
                            1, soda_reviews)
    >>> restaurant_name(soda)
    'Soda'
    >>> restaurant_location(soda)
    [127.0, 0.1]
    >>> restaurant_ratings(soda)
    [4.5, 4]
    '''
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    return restaurant["Name"]
    # END Question 2

def restaurant_location(restaurant):
    """Return the location of the restaurant, which is a list containing
    latitude and longitude."""
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    return restaurant["Location"]
    # END Question 2

def restaurant_categories(restaurant):
    """Return the categories of the restaurant, which is a list of strings."""
    '''
    Q: Given a restaurant, what does restaurant_categories in
    abstractions.py return?
    Choose the number of the correct choice:
    0) a single number (rating)
    1) a list of numbers (ratings)
    2) a single string (category)
    3) a list of strings (categories)
    3

    '''
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    return restaurant["Categories"]
    # END Question 2

def restaurant_price(restaurant):
    """Return the price of the restaurant, which is a number."""
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    return restaurant["Price"]
    # END Question 2

def restaurant_ratings(restaurant):
    """Return a list of ratings, which are numbers from 1 to 5, of the
    restaurant based on the reviews of the restaurant."""
    # BEGIN Question 2
    "*** REPLACE THIS LINE ***"
    return restaurant["Reviews"]
    # END Question 2