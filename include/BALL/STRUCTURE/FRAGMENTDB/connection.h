// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_CONNECTION_H
#define BALL_STRUCTURE_FRAGMENTDB_CONNECTION_H

#include<BALL/KERNEL/atom.h>
#include<BALL/KERNEL/bond.h>
#include<BALL/DATATYPE/string.h>
#include<BALL/CONCEPT/persistentObject.h>

namespace BALL {

	struct Connection : public PersistentObject
	{
		Atom*				atom;
		String			type_name;
		String			connect_to;
		Bond::Order order;
		float				dist;
		float				delta;
		Connection& operator= (const Connection& other)
		{
			atom = other.atom;
			type_name = other.type_name;
			connect_to = other.connect_to;
			order = other.order;
			dist = other.dist;
			delta = other.delta;
			return *this;
		}
		Connection(const Connection& other)
			: PersistentObject(other),
			atom(other.atom),
			type_name ( other.type_name),
			connect_to ( other.connect_to),
			order ( other.order),
			dist ( other.dist),
			delta ( other.delta)
		{
		}
		Connection()
			: atom(0),
			type_name (""),
			connect_to (""),
			order (0),
			dist (0.f),
			delta (0.f)
		{
		}
		~Connection()
		{
		asm("nop");
		}
	};

}

#endif // CONNECTION_H
