import sqlite3
import os.path
import pandas as pd
import numpy as np


def create_db():
	if os.path.exists("players.db"):
		os.remove("players.db")
	data_file = 'players.csv'
	fifa_players = pd.read_csv(data_file)

	c = sqlite3.connect('players.db')

	# Create and populate the Bio table
	c = sqlite3.connect('players.db')
	c.execute('''CREATE TABLE Bio(
		ID numeric PRIMARY KEY, full_name text,
		club text, nationality text,
		birth_date text, age numeric,
		height_cm numeric, weight_kg numeric,
		photo text)''')

	subset = fifa_players[['ID', 'full_name', 'club', 'nationality', 'birth_date', 'age',
		'height_cm', 'weight_kg', 'photo']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Bio (ID, full_name, club,nationality, birth_date, age,
			height_cm, weight_kg, photo) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the the Overall_stat relation
	c.execute('''CREATE TABLE Overall_stat(
		playerID numeric, preferred_foot text,
		overall numeric, potential numeric,
		pac numeric, body_type text,
		weak_foot numeric, international_reputation numeric,
		stamina numeric, strength numeric,
		balance numeric, reactions numeric,
		heading_accuracy numeric, interceptions numeric,
		positioning numeric, vision numeric,
		penalties numeric, composure numeric,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'preferred_foot', 'overall', 'potential', 'pac', 'body_type',
		'weak_foot', 'international_reputation', 'stamina', 'strength', 'balance', 'reactions', 'heading_accuracy',
		'interceptions', 'positioning', 'vision', 'penalties', 'composure']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Overall_stat (playerID, preferred_foot, overall, potential, pac, body_type,
			weak_foot, international_reputation, stamina, strength, balance, reactions, heading_accuracy,
			interceptions, positioning, vision, penalties, composure)
			VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))


	# Create and populate the Attacking_stat relation
	c.execute('''CREATE TABLE Attacking_stat(
		playerID numeric, sho numeric,
		pas numeric, dri numeric,
		crossing numeric, finishing numeric,
		short_passing numeric, volleys numeric,
		dribbling numeric, curve numeric,
		free_kick_accuracy numeric, long_passing numeric,
		ball_control numeric, acceleration numeric,
		sprint_speed numeric, agility numeric,
		shot_power numeric, long_shots numeric,
		aggression numeric, work_rate_att text,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'sho', 'pas', 'dri', 'crossing', 'finishing', 'short_passing',
	'volleys', 'dribbling', 'curve', 'free_kick_accuracy', 'long_passing', 'ball_control',
	'acceleration', 'sprint_speed', 'agility', 'shot_power', 'long_shots', 'aggression', 'work_rate_att']]

	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Attacking_stat (playerID, sho, pas, dri, crossing, finishing, short_passing,
			volleys, dribbling, curve, free_kick_accuracy, long_passing, ball_control,
			acceleration, sprint_speed, agility, shot_power, long_shots, aggression, work_rate_att)
			VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Defense_stat relation
	c.execute('''CREATE TABLE Defense_stat(
		playerID numeric, def numeric,
		phy numeric, work_rate_def text,
		marking numeric, standing_tackle numeric,
		sliding_tackle numeric,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'def', 'phy', 'work_rate_def', 'marking', 'standing_tackle', 'sliding_tackle']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Defense_stat (playerID, def, phy, work_rate_def, marking, standing_tackle, sliding_tackle)
			VALUES (?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Goalkeeper_stat relation
	c.execute('''CREATE TABLE Goalkeeper_stat(
		playerID numeric, gk_diving numeric,
		gk_handling numeric, gk_kicking text,
		gk_positioning numeric, gk_reflexes numeric,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'gk_diving', 'gk_handling', 'gk_kicking', 'gk_positioning', 'gk_reflexes']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Goalkeeper_stat (playerID, gk_diving, gk_handling, gk_kicking, gk_positioning, gk_reflexes)
			VALUES (?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Money relation
	c.execute('''CREATE TABLE Money(
		playerID numeric, eur_wage numeric,
		eur_value numeric, eur_release_clause text,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'eur_wage', 'eur_value', 'eur_release_clause']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Money (playerID, eur_wage, eur_value, eur_release_clause)
			VALUES (?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Teams relation
	c.execute('''CREATE TABLE Teams(
		playerID numeric, club text,
		league text, clublogo text,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'club', 'league', 'club_logo']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Teams (playerID, club, league, clublogo)
			VALUES (?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Nationality relation
	c.execute('''CREATE TABLE Nationality(
		playerID numeric, country_name text,
		flag text, FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'nationality', 'flag']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Nationality (playerID, country_name, flag)
			VALUES (?, ?, ?)""", list(row_cast))

	# Create and populate the Positions relation
	c.execute('''CREATE TABLE Positions(
		playerID numeric, rs numeric,
		rw numeric, rf numeric,
		ram numeric, rdm numeric,
		rcb numeric, rm numeric,
		rb numeric, rwb numeric,
		cf numeric, cam numeric,
		cdm numeric, cm numeric,
		cb numeric, ls numeric,
		lw numeric, lf numeric,
		lam numeric, ldm numeric,
		lcb numeric, lm numeric,
		lb numeric, lwb numeric,
		gk numeric, FOREIGN KEY(playerID) REFERENCES Bio(ID))''')

	subset = fifa_players[['ID', 'rs', 'rw', 'rf', 'ram', 'rdm', 'rcb', 'rm', 'rb', 'rwb',
		'cf', 'cam', 'cdm', 'cm', 'cb', 'ls', 'lw', 'lf', 'lam', 'ldm', 'lcb', 'lm', 'lb',
		'lwb', 'gk']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		c.execute("""INSERT INTO Positions (playerID, rs, rw, rf, ram, rdm, rcb, rm, rb, rwb,
		cf, cam, cdm, cm, cb, ls, lw, lf, lam, ldm, lcb, lm, lb, lwb, gk)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Trait relation
	c.execute('''CREATE TABLE Trait(
		playerID numeric, chip_shot_trait numeric,
		corner_specialist_trait numeric, diver_trait numeric,
		finesse_shot_trait numeric, gk_long_throw_trait numeric,
		gk_up_for_corners_trait numeric, injury_free_trait numeric,
		injury_prone_trait numeric, leadership_trait numeric,
		long_passer_trait numeric, long_shot_taker_trait numeric,
		one_club_player_trait numeric, playmaker_trait numeric,
		power_free_kick_trait numeric, power_header_trait numeric,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')


	subset = fifa_players[['ID', 'chip_shot_trait', 'corner_specialist_trait', 'diver_trait',
		'finesse_shot_trait', 'gk_long_throw_trait', 'gk_up_for_corners_trait', 'injury_free_trait',
		'injury_prone_trait', 'leadership_trait', 'long_passer_trait', 'long_shot_taker_trait',
		'one_club_player_trait', 'playmaker_trait', 'power_free_kick_trait', 'power_header_trait']]
	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		row_cast = [x if (type(x) is not np.bool_) else bool(x) for x in row_cast]
		if i == 0:
			print(row_cast)
		c.execute("""INSERT INTO Trait (playerID, chip_shot_trait, corner_specialist_trait, diver_trait,
		finesse_shot_trait, gk_long_throw_trait, gk_up_for_corners_trait, injury_free_trait,
		injury_prone_trait, leadership_trait, long_passer_trait, long_shot_taker_trait,
		one_club_player_trait, playmaker_trait, power_free_kick_trait, power_header_trait)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	# Create and populate the Specialty relation
	c.execute('''CREATE TABLE Specialty(
		playerID numeric, speedster_speciality numeric,
		dribbler_speciality numeric, engine_speciality numeric,
		distance_shooter_speciality numeric, free_kick_specialist_speciality numeric,
		tackling_speciality numeric, strength_speciality numeric,
		FOREIGN KEY(playerID) REFERENCES Bio(ID))''')


	subset = fifa_players[['ID', 'speedster_speciality', 'dribbler_speciality', 'engine_speciality',
		'distance_shooter_speciality', 'free_kick_specialist_speciality', 'tackling_speciality',
		'strength_speciality']]

	for i in subset.axes[0]:
		row_cast = [x if (type(x) is not np.int64) else int(x) for x in subset.iloc[i,:]]
		row_cast = [x if (type(x) is not np.bool_) else bool(x) for x in row_cast]
		c.execute("""INSERT INTO Specialty (playerID, speedster_speciality, dribbler_speciality, engine_speciality,
		distance_shooter_speciality, free_kick_specialist_speciality, tackling_speciality,
		strength_speciality)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?)""", list(row_cast))

	c.commit()
	c.close()
	return "OK"

if __name__ == "__main__":
	create_db()