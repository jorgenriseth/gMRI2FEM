import simple_mri as sm
import numpy as np
import matplotlib.pyplot as plt
import skimage


def slice_volume(volume, orientation, idx):
    tail_dims = list(range(3, volume.ndim))
    if orientation == "sagittal":
        slice_ = np.s_[idx]
    elif orientation == "coronal":
        slice_ = np.s_[:, idx]
    elif orientation in ["axial", "transversal"]:
        slice_ = np.s_[:, :, idx]
    else:
        raise ValueError(f"invalid orientation {orientation}")
    return np.rot90(volume[slice_])


def crop_limits(im, crop_margin):
    try:
        thresh = skimage.filters.threshold_triangle(np.where(np.isnan(im), 0, im))
    except ValueError:
        return 0, im.shape[0], 0, im.shape[1]
    binary = im > thresh
    binary = skimage.filters.gaussian(binary, 5)
    try:
        binary = binary > skimage.filters.threshold_isodata(binary)
    except IndexError:
        return 0, im.shape[0], 0, im.shape[1]
    bounds_i, bounds_j = [(min(x), max(x)) for x in np.where(binary)]
    crop_i = int(max(0, bounds_i[0] - crop_margin)), int(bounds_i[1] + crop_margin)
    crop_j = int(max(0, bounds_j[0] - crop_margin)), int(bounds_j[1] + crop_margin)
    return *crop_i, *crop_j


def compute_row_fractional_widths(
    row: list[np.ndarray],
) -> dict[str, tuple[float, ...]]:
    heights, widths = zip(*[im.shape for im in row])
    H = heights[0]
    scaled_widths = tuple((H / hi) * wi for hi, wi in zip(heights, widths))
    W = sum(scaled_widths)
    fractional_widths = tuple(round(swi / W, 3) for swi in scaled_widths)
    return {"shape": (H, W), "fractional_column_widths": fractional_widths}


def compute_body_fractional_heights(images: list[list[np.ndarray]]):
    row_sizes = [compute_row_fractional_widths(row) for row in images]
    heights, widths = zip(*[row_size["shape"] for row_size in row_sizes])
    W = widths[0]
    scaled_heights = tuple((W / wi) * hi for hi, wi in zip(heights, widths))

    H = sum(scaled_heights)
    fractional_widths = tuple(row["fractional_column_widths"] for row in row_sizes)
    fractional_heights = tuple(round(shi / H, 3) for shi in scaled_heights)
    horizontal_offsets = [
        (0, *np.cumsum(row["fractional_column_widths"])) for row in row_sizes
    ]
    vertical_offsets = [round(1 - hi, 3) for hi in np.cumsum(fractional_heights)]
    body_axis_sizes = [
        [(wj, hi) for wj in row["fractional_column_widths"]]
        for hi, row in zip(fractional_heights, row_sizes)
    ]
    body_axis_corners = [
        [(xj, yi) for xj in row_horizontal_offsets]
        for yi, row_horizontal_offsets in zip(vertical_offsets, horizontal_offsets)
    ]

    return {"shape": (H, W), "offsets": body_axis_corners, "sizes": body_axis_sizes}


def build_image_grid(images, grid_config, frame=False):
    fig_width = grid_config["fig_width"]
    body_dims = compute_body_fractional_heights(images)
    sizes, offsets = body_dims["sizes"], body_dims["offsets"]
    body_height, body_width = body_dims["shape"]
    body_aspect = body_width / body_height

    left_margin = grid_config["fig_left_margin"] + grid_config["row_left_margin"]
    right_margin = grid_config["fig_right_margin"] + grid_config["row_right_margin"]
    top_margin = grid_config["fig_top_margin"] + grid_config["col_top_margin"]
    bottom_margin = grid_config["fig_bottom_margin"] + grid_config["col_bottom_margin"]

    body_width = 1 - sum(
        grid_config[key]
        for key in (
            "fig_left_margin",
            "row_left_margin",
            "row_right_margin",
            "fig_right_margin",
        )
    )
    body_height = 1 - sum(
        grid_config[key]
        for key in (
            "fig_top_margin",
            "col_top_margin",
            "col_bottom_margin",
            "fig_bottom_margin",
        )
    )

    left_column_sizes = [row_sizes[0] for row_sizes in sizes]
    left_column_offsets = [row_offsets[0] for row_offsets in offsets]
    right_column_sizes = [row_sizes[-1] for row_sizes in sizes]
    right_column_offsets = [row_offsets[-1] for row_offsets in offsets]

    fig_aspect = body_aspect / (body_width / body_height)
    fig_height = fig_width / fig_aspect
    fig = plt.figure(figsize=(fig_width, fig_height))
    fig_margins = {
        "fig_top": fig.add_axes(
            (
                left_margin,
                1 - grid_config["fig_top_margin"],
                body_width,
                grid_config["fig_top_margin"],
            ),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_bottom": fig.add_axes(
            (left_margin, 0, body_width, grid_config["fig_bottom_margin"]),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_left": fig.add_axes(
            (0, bottom_margin, grid_config["fig_left_margin"], body_height),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_right": fig.add_axes(
            (
                1 - grid_config["fig_right_margin"],
                bottom_margin,
                grid_config["fig_right_margin"],
                body_height,
            ),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
    }
    inner_margins = {
        "col_top": [
            fig.add_axes(
                (
                    left_margin + body_width * oi,
                    1 - top_margin,
                    body_width * wi,
                    grid_config["col_top_margin"],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for (oi, _), (wi, _) in zip(offsets[0], sizes[0])
        ],
        "col_bottom": [
            fig.add_axes(
                (
                    left_margin + body_width * oi,
                    grid_config["fig_bottom_margin"],
                    body_width * wi,
                    grid_config["col_bottom_margin"],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for (oi, _), (wi, _) in zip(offsets[-1], sizes[-1])
        ],
        "row_left": [
            fig.add_axes(
                (
                    grid_config["fig_left_margin"],
                    bottom_margin + body_height * oi,
                    grid_config["row_left_margin"],
                    body_height * hi,
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for (_, oi), (wi, hi) in zip(left_column_offsets, left_column_sizes)
        ],
        "row_right": [
            fig.add_axes(
                (
                    1 - right_margin,
                    bottom_margin + body_height * oi,
                    grid_config["row_right_margin"],
                    body_height * hi,
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for (_, oi), (wi, hi) in zip(right_column_offsets, right_column_sizes)
        ],
    }
    body_axes = {
        "body": [
            [
                fig.add_axes(
                    (
                        left_margin + body_width * oi,
                        bottom_margin + body_height * vi,
                        body_width * wi,
                        body_height * hi,
                    ),
                    xticks=[],
                    yticks=[],
                    frameon=True,
                )
                for (oi, vi), (wi, hi) in zip(row_offsets, row_sizes)
            ]
            for row_offsets, row_sizes in zip(offsets, sizes)
        ]
    }

    axes = {
        **fig_margins,
        **inner_margins,
        **body_axes,
    }
    return fig, axes


def legacy_image_grid_builder(fig_width, im_width, im_height, grid_config, frame=False):
    body_width = grid_config["ncols"] * im_width
    body_height = grid_config["nrows"] * im_height

    col_widths = [
        *[
            grid_config[prop] * body_width
            for prop in ["fig_left_margin", "row_left_margin"]
        ],
        *([im_width] * grid_config["ncols"]),
        *[
            grid_config[prop] * body_width
            for prop in ["row_right_margin", "fig_right_margin"]
        ],
    ]
    row_heights = [
        *[
            grid_config[prop] * body_height
            for prop in ["fig_top_margin", "col_top_margin"]
        ],
        *([im_height] * grid_config["nrows"]),
        *[
            grid_config[prop] * body_height
            for prop in ["col_bottom_margin", "fig_bottom_margin"]
        ],
    ]

    total_width = sum(col_widths)
    total_height = sum(row_heights)
    aspect_ratio = total_width / total_height
    fig_height = fig_width / aspect_ratio

    col_widths = np.array(col_widths) / total_width
    row_heights = np.array(row_heights) / total_height
    col_offsets = np.cumsum([0] + list(col_widths[:-1]))
    row_offsets = 1 - np.cumsum(row_heights)

    fig = plt.figure(figsize=(fig_width, fig_height))
    axes = {
        "fig_top": fig.add_axes(
            (col_offsets[2], row_offsets[0], sum(col_widths[2:-2]), row_heights[0]),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_bottom": fig.add_axes(
            (col_offsets[2], row_offsets[-1], sum(col_widths[2:-2]), row_heights[-1]),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_left": fig.add_axes(
            (col_offsets[0], row_offsets[-3], col_widths[0], sum(row_heights[2:-2])),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "fig_right": fig.add_axes(
            (col_offsets[-1], row_offsets[-3], col_widths[-1], sum(row_heights[2:-2])),
            xticks=[],
            yticks=[],
            frameon=frame,
        ),
        "col_top": [
            fig.add_axes(
                (
                    col_offsets[2 + idx],
                    row_offsets[1],
                    col_widths[2 + idx],
                    row_heights[1],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for idx in range(grid_config["ncols"])
        ],
        "col_bottom": [
            fig.add_axes(
                (
                    col_offsets[2 + idx],
                    row_offsets[-2],
                    col_widths[2 + idx],
                    row_heights[-2],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for idx in range(grid_config["ncols"])
        ],
        "row_left": [
            fig.add_axes(
                (
                    col_offsets[1],
                    row_offsets[2 + idx],
                    col_widths[1],
                    row_heights[2 + idx],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for idx in range(grid_config["nrows"])
        ],
        "row_right": [
            fig.add_axes(
                (
                    col_offsets[-2],
                    row_offsets[2 + idx],
                    col_widths[-2],
                    row_heights[2 + idx],
                ),
                xticks=[],
                yticks=[],
                frameon=frame,
            )
            for idx in range(grid_config["nrows"])
        ],
        "body": [
            fig.add_axes(
                (
                    col_offsets[2 + col],
                    row_offsets[2 + row],
                    col_widths[2 + col],
                    row_heights[2 + row],
                ),
                xticks=[],
                yticks=[],
            )
            for row, col in map(
                lambda idx: (idx // grid_config["ncols"], idx % grid_config["ncols"]),
                range(grid_config["nrows"] * grid_config["ncols"]),
            )
        ],
    }
    return fig, axes


def flatten(list_of_lists):
    return sum(list_of_lists, start=[])
