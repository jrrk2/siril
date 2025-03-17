/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/siril.h"
#include "gui/user_polygons.h"

#define FROM_BE64_INTO(dest, val, type) \
do { \
	union { type v; uint64_t i; } conv; \
	memcpy(&conv.i, &val, sizeof(type)); \
	conv.i = GUINT64_FROM_BE(conv.i); \
	(dest) = conv.v; \
} while(0)

int get_unused_polygon_id(void) {
	int id = 1;
	GList *l;

	while (TRUE) {
		gboolean id_in_use = FALSE;
		for (l = gui.user_polygons; l != NULL; l = l->next) {
			UserPolygon *polygon = (UserPolygon *)l->data;
			if (polygon->id == id) {
				id_in_use = TRUE;
				break;
			}
		}

		if (!id_in_use) {
			return id;
		}

		id++;
	}
}

int add_user_polygon(point *points, int n_points, const GdkRGBA *color, gboolean fill) {
	int id = get_unused_polygon_id();

	UserPolygon *polygon = g_new(UserPolygon, 1);
	polygon->id = id;
	polygon->n_points = n_points;
	polygon->points = g_new(point, n_points);
	polygon->fill = fill;
	polygon->color = *color;

	for (int i = 0; i < n_points; i++) {
		polygon->points[i] = points[i];
	}

	gui.user_polygons = g_list_append(gui.user_polygons, polygon);

	return id; // Return the generated ID to the caller
}

static void free_user_polygon(gpointer data) {
	UserPolygon *polygon = (UserPolygon *)data;
	if (polygon) {
		g_free(polygon->points);
		g_free(polygon);
	}
}

gboolean delete_user_polygon(int id) {
	GList *l;
	for (l = gui.user_polygons; l != NULL; l = l->next) {
		UserPolygon *polygon = (UserPolygon *)l->data;
		if (polygon->id == id) {
			gui.user_polygons = g_list_delete_link(gui.user_polygons, l);
			free_user_polygon(polygon);
			return TRUE;
		}
	}
	return FALSE;
}

void clear_user_polygons(void) {
	g_list_free_full(gui.user_polygons, free_user_polygon);
	gui.user_polygons = NULL;
}

UserPolygon* deserialize_polygon(const uint8_t *data, size_t size) {
	if (size < 13) {  // 4 + 4 + 4 + 1 bytes (minimum)
		fprintf(stderr, "Invalid data size\n");
		return NULL;
	}

	const uint8_t *ptr = data;
	UserPolygon *polygon = g_malloc0(sizeof(UserPolygon));
	if (!polygon) {
		fprintf(stderr, "Memory allocation failed for polygon\n");
		return NULL;
	}

	// Read ID
	polygon->id = (int32_t)g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	// Read number of points
	polygon->n_points = (int32_t)g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	// Read RGBA (packed into uint32_t)
	uint32_t packed_color = g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	polygon->color.red   = ((packed_color >> 24) & 0xFF) / 255.0;
	polygon->color.green = ((packed_color >> 16) & 0xFF) / 255.0;
	polygon->color.blue  = ((packed_color >> 8)  & 0xFF) / 255.0;
	polygon->color.alpha = (packed_color & 0xFF) / 255.0;

	// Read fill flag (1 byte boolean)
	polygon->fill = (gboolean)(*ptr);
	ptr += 1;  // Move pointer past boolean

	if (polygon->n_points <= 0 || polygon->n_points > MAX_POLYGON_POINTS) {
		fprintf(stderr, "Invalid number of points: %d\n", polygon->n_points);
		g_free(polygon);
		return NULL;
	}

	// Allocate memory for points
	polygon->points = g_malloc0(polygon->n_points * sizeof(point));
	if (!polygon->points) {
		fprintf(stderr, "Memory allocation for points failed\n");
		g_free(polygon);
		return NULL;
	}

	// Ensure enough data for points
	size_t required_size = 13 + polygon->n_points * (2 * sizeof(double));
	if (size < required_size) {
		fprintf(stderr, "Not enough data for points\n");
		g_free(polygon->points);
		g_free(polygon);
		return NULL;
	}

	// Read points
	for (int i = 0; i < polygon->n_points; i++) {
		double x_BE, y_BE;
		memcpy(&x_BE, ptr, sizeof(double));
		ptr += sizeof(double);
		FROM_BE64_INTO(polygon->points[i].x, x_BE, double);

		memcpy(&y_BE, ptr, sizeof(double));
		ptr += sizeof(double);
		FROM_BE64_INTO(polygon->points[i].y, y_BE, double);
	}

	return polygon;
}
