/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2019
 */

#pragma once


#include <memory>

#if __cplusplus >= 201402L

#define MEXOCT_BUFFER(T, buf, size)                               \
  auto __mexoct_buffer_ ## buf = std::make_unique<T []> (size);     \
  T *buf = __mexoct_buffer__ ## buf.get ()

#else

#define MEXOCT_BUFFER(T, buf, size)                               \
  std::unique_ptr<T []> __mexoct_buffer__ ## buf { new T [size] };   \
  T *buf = __mexoct_buffer__ ## buf.get ()
#endif


